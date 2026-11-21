// Prevent console window in addition to Slint window in Windows release builds when, e.g., starting the app via file manager. Ignored on other platforms.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use femodel::prelude::*;
use femodel::utils::generations::{
    FieldType, Perturbation, RandomField, save_fields_to_abaqus, save_fields_to_lsdyna,
};
use rfd::AsyncFileDialog;
use slint::{Model, ModelRc, StandardListViewItem, VecModel};
use std::cell::RefCell;
use std::error::Error;
use std::rc::Rc;
use std::sync::{Arc, RwLock};

#[derive(Clone)]
struct FieldDraft {
    id: usize,
    name: String,
    mean: f32,
    std: f32,
    inverse: bool,
}

#[derive(Clone)]
struct PerturbationDraft {
    field_type: FieldType,
    fields: Vec<FieldDraft>,
}

struct AppState {
    model: Arc<RwLock<FEModel>>, // None пока файл не открыт
    perturbations: Vec<PerturbationDraft>,
}

/// Всё, что нужно фоновому потоку. Только `Send`-данные — ни `Rc`, ни `AppState`.
struct GenerationJob {
    model: Arc<RwLock<FEModel>>,
    drafts: Vec<PerturbationDraft>,
    output: std::path::PathBuf,
    save_vtk: bool,
}

impl AppState {
    /// Клонирует данные для фонового потока. Вызывать только на UI-потоке.
    fn snapshot_for_generation(&self, output: std::path::PathBuf, save_vtk: bool) -> GenerationJob {
        GenerationJob {
            model: Arc::clone(&self.model),
            drafts: self.perturbations.clone(),
            output,
            save_vtk,
        }
    }

    fn ready_to_generate(&self) -> bool {
        let has_nodes = self
            .model
            .read()
            .map(|m| !m.nodes.is_empty())
            .unwrap_or(false);
        has_nodes && self.perturbations.iter().any(|p| !p.fields.is_empty())
    }
}

/// Генерирует и записывает файлы. Чистая работа, выполняется в фоновом потоке.
/// Возвращает текст статуса для UI.
fn run_generation(job: GenerationJob) -> String {
    let model = job.model.read().unwrap();

    let results: Vec<RandomField> = job
        .drafts
        .iter()
        .flat_map(|p| {
            let pert = Perturbation::new(&model, p.field_type);
            p.fields.iter().map(move |f| {
                RandomField::new()
                    .with_id(f.id)
                    .with_name(f.name.clone())
                    .generate(&pert, f.mean, f.std, f.inverse)
            })
        })
        .collect();

    let refs: Vec<&RandomField> = results.iter().collect();
    let path = job.output.to_string_lossy();
    let ext = job
        .output
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    let mut status = match ext.as_str() {
        "k" => save_fields_to_lsdyna(&path, &refs),
        _ => save_fields_to_abaqus(&path, &refs),
    }
    .map(|()| format!("Сохранено: {}", job.output.display()))
    .unwrap_or_else(|e| format!("Ошибка записи {}: {e}", job.output.display()));

    if job.save_vtk {
        let vtk_path = job.output.with_extension("vtk");
        let names: Vec<String> = results.iter().map(|f| f.name.clone()).collect();
        let nodal_fields: Vec<NodeScalarField> = names
            .iter()
            .zip(results.iter())
            .map(|(name, f)| (name.as_str(), &f.data))
            .collect();
        match model.save_vtk(
            vtk_path.to_string_lossy().as_ref(),
            Some(&nodal_fields),
            None,
        ) {
            Ok(()) => status += &format!("\nVTK: {}", vtk_path.display()),
            Err(e) => status += &format!("\nОшибка VTK: {e}"),
        }
    }
    status
}

slint::include_modules!();

fn main() -> Result<(), Box<dyn Error>> {
    let ui = AppWindow::new()?;
    let ui_handle = ui.as_weak();
    let app_state = Rc::new(RefCell::new(AppState {
        model: Arc::new(RwLock::new(FEModel::default())),
        perturbations: Vec::new(),
    }));

    // let m = Rc::new(std::cell::RefCell::new(FEModel::default()));

    let items: Vec<StandardListViewItem> = Vec::new();
    let perturbation_model = Rc::new(VecModel::from(items));
    ui.set_perturbations(ModelRc::from(perturbation_model.clone()));

    let items: Vec<StandardListViewItem> = Vec::new();
    let fields_model = Rc::new(VecModel::from(items));
    ui.set_random_fields(ModelRc::from(fields_model.clone()));

    // выбор и чтение файла с сеточной моделью
    // let app = app_state.clone();
    // let ui_handle1 = ui_handle.clone();
    ui.on_choose_file({
        let app_state = Rc::clone(&app_state);
        let ui_handle = ui_handle.clone();
        move || {
            let app_state = app_state.clone();
            let ui_handle = ui_handle.clone();
            slint::spawn_local(async move {
                let pick = AsyncFileDialog::new()
                    .add_filter("supported meshes", &["k", "inp"])
                    .add_filter("ls-dyna mesh", &["k"])
                    .add_filter("abaqus mesh", &["inp"])
                    .pick_file()
                    .await;
                if let Some(path) = pick {
                    let path = String::from(path.path().to_string_lossy());
                    println!("Reading model...");
                    app_state
                        .borrow_mut()
                        .model
                        .write()
                        .unwrap()
                        .read(&path)
                        .unwrap(); // Сохраняем в UI
                    let st = app_state.borrow();
                    let m = st.model.read().unwrap();
                    println!("{}", m);
                    let ncount = m.nodes.len();
                    let solids = m.solids.len();
                    let shells = m.shells.len();
                    let bb = m.bbox();
                    drop(m);
                    drop(st);
                    let mut s = format!(
                        "Узлов: {ncount}, объемных элементов: {solids}, оболочек: {shells}\n"
                    );
                    s += &format!(
                        "Габариты модели: x: {}->{}, y: {}->{}, z: {}->{}",
                        bb.min_x, bb.max_x, bb.min_y, bb.max_y, bb.min_z, bb.max_z
                    );
                    let _ = ui_handle.upgrade_in_event_loop(move |ui| {
                        ui.set_file_path(path.into());
                        ui.set_model_info(s.into());
                    });
                }
            })
            .unwrap();
        }
    });

    ui.on_add_perturbation({
        let app_state = Rc::clone(&app_state);
        let ui = ui_handle.unwrap();
        let pm = Rc::clone(&perturbation_model);
        move |field_type: i32, params: ModelRc<f32>| {
            let params: Vec<f32> = params.iter().collect();
            // let perturbation_model = pm.clone();
            match field_type {
                0 => {
                    let lc = params.first().unwrap_or(&1.0);
                    let p = PerturbationDraft {
                        field_type: FieldType::Isotropic(*lc),
                        fields: Vec::new(),
                    };
                    app_state.borrow_mut().perturbations.push(p);
                    pm.push(StandardListViewItem::from(
                        String::from(&FieldType::Isotropic(*lc)).as_str(),
                    ));
                }
                1 => {
                    let p = PerturbationDraft {
                        field_type: FieldType::Orthotropic(params[0], params[1], params[2]),
                        fields: Vec::new(),
                    };
                    app_state.borrow_mut().perturbations.push(p);
                    pm.push(StandardListViewItem::from(
                        String::from(&FieldType::Orthotropic(params[0], params[1], params[2]))
                            .as_str(),
                    ));
                }
                _ => {}
            }
            let idx = (app_state.borrow().perturbations.len() - 1) as i32;
            ui.set_current_perturbation_index(idx);
        }
    });

    // let app = app_state.clone();
    // let pm = perturbation_model.clone();
    ui.on_delete_current_perturbation({
        let app_state = Rc::clone(&app_state);
        let pm = Rc::clone(&perturbation_model);
        move |idx: i32| {
            let n = app_state.borrow().perturbations.len();
            if idx == -1 || n == 0 || idx as usize >= n {
                return;
            };
            // let perturbation_model = pm;
            app_state.borrow_mut().perturbations.remove(idx as usize);
            pm.remove(idx as usize);
        }
    });

    ui.on_add_field({
        let app_state = app_state.clone();
        let fm = fields_model.clone();
        move |idx, fp: FieldProps| {
            app_state.borrow_mut().perturbations[idx as usize]
                .fields
                .push(FieldDraft {
                    id: fp.id as usize,
                    name: String::from(&fp.name),
                    mean: fp.mean,
                    std: fp.std,
                    inverse: fp.inverse,
                });
            fm.push(StandardListViewItem::from(
                String::from(format!(
                    "id({}), name({}), mean({}), std({}), inverse({})",
                    fp.id, fp.name, fp.mean, fp.std, fp.inverse
                ))
                .as_str(),
            ));
        }
    });

    ui.on_change_current_perturbation({
        let app_state = app_state.clone();
        let fm = fields_model.clone();
        move |idx| {
            fm.clear();
            let n = app_state.borrow().perturbations.len();
            if idx == -1 || n == 0 || idx as usize >= n {
                return;
            }
            for f in app_state.borrow().perturbations[idx as usize].fields.iter() {
                fm.push(StandardListViewItem::from(
                    String::from(format!(
                        "id({}), name({}), mean({}), std({}), inverse({})",
                        f.id, f.name, f.mean, f.std, f.inverse
                    ))
                    .as_str(),
                ));
            }
        }
    });

    ui.on_delete_current_field({
        let app_state = app_state.clone();
        let fm = fields_model.clone();
        move |p_idx, f_idx| {
            app_state.borrow_mut().perturbations[p_idx as usize]
                .fields
                .remove(f_idx as usize);
            fm.remove(f_idx as usize);
        }
    });

    ui.on_generate({
        let app_state = app_state.clone();
        let ui_handle = ui_handle.clone();
        let ui = ui_handle.unwrap();

        move || {
            let snapshot = {
                let st = app_state.borrow();
                if !st.ready_to_generate() {
                    return;
                }
                st
            };

            // Читаем UI-свойство и готовим снимок — оба синхронно, на UI-потоке.
            let save_vtk = ui.get_is_save_vtk();
            let model_and_drafts = (Arc::clone(&snapshot.model), snapshot.perturbations.clone());
            drop(snapshot);

            ui_handle
                .upgrade_in_event_loop(|ui| ui.set_generating(true))
                .ok();

            let ui_handle = ui_handle.clone();
            // log::info!("generation started");
            // let app_state = app_state.clone();

            slint::spawn_local(async move {
                let Some(file) = AsyncFileDialog::new()
                    .add_filter("ls-dyna", &["k"])
                    .add_filter("abaqus", &["inp"])
                    .set_file_name("fields")
                    .save_file()
                    .await
                else {
                    ui_handle
                        .upgrade_in_event_loop(|ui| ui.set_generating(false))
                        .ok();
                    return;
                };

                let job = {
                    // let st = app_state.borrow();
                    let (model, drafts) = &model_and_drafts;
                    GenerationJob {
                        model: Arc::clone(model),
                        drafts: drafts.clone(),
                        output: file.path().to_path_buf(),
                        save_vtk,
                    }
                };

                let ui_handle = ui_handle.clone();
                std::thread::spawn(move || {
                    let status = run_generation(job);
                    ui_handle
                        .upgrade_in_event_loop(move |ui| {
                            ui.set_model_info(status.into());
                            ui.set_generating(false);
                        })
                        .ok();
                });
            })
            .ok();
        }
    });

    ui.on_get_current_field_data({
        // let ui = ui_handle.unwrap();
        let app_state = app_state.clone();
        move |p_idx, f_idx| {
            let app_state = app_state.borrow();
            if let Some(p) = app_state.perturbations.get(p_idx as usize) {
                if let Some(f) = p.fields.get(f_idx as usize) {
                    return FieldProps {
                        id: f.id as i32,
                        name: String::from(&f.name).into(),
                        mean: f.mean,
                        std: f.std,
                        inverse: f.inverse,
                    };
                }
            };
            FieldProps {
                id: 1,
                name: "".into(),
                mean: 0.0,
                std: 1.0,
                inverse: false,
            }
        }
    });

    ui.on_edit_field({
        let app_state = app_state.clone();
        let fm = fields_model.clone();
        move |p_idx, f_idx, fp| {
            let mut app_state = app_state.borrow_mut();
            if let Some(p) = app_state.perturbations.get_mut(p_idx as usize) {
                if let Some(f) = p.fields.get_mut(f_idx as usize) {
                    f.id = fp.id as usize;
                    f.name = String::from(&fp.name);
                    f.mean = fp.mean;
                    f.std = fp.std;
                    f.inverse = fp.inverse;
                }
            };
            let item = StandardListViewItem::from(
                format!(
                    "id({}), name({}), mean({}), std({}), inverse({})",
                    fp.id, fp.name, fp.mean, fp.std, fp.inverse
                )
                .as_str(),
            );
            fm.set_row_data(f_idx as usize, item);
        }
    });

    ui.on_get_current_perturbation_data({
        let app_state = app_state.clone();
        move |p_idx: i32| -> PerturbationProps {
            let (kind, params): (i32, Vec<f32>) = {
                let st = app_state.borrow();
                match st.perturbations.get(p_idx as usize).map(|p| p.field_type) {
                    Some(FieldType::Isotropic(lc)) => (0, vec![lc]),
                    Some(FieldType::Orthotropic(lx, ly, lz)) => (1, vec![lx, ly, lz]),
                    None => (0, vec![1.0]),
                }
            };
            PerturbationProps {
                kind,
                params: ModelRc::from(Rc::new(VecModel::from(params))),
            }
        }
    });

    ui.on_edit_perturbation({
        let app_state = app_state.clone();
        let pm = perturbation_model.clone();
        move |idx, pprops| {
            let mut st = app_state.borrow_mut();
            if let Some(p) = st.perturbations.get_mut(idx as usize) {
                let params: Vec<f32> = pprops.params.iter().collect();
                let s: String = match pprops.kind {
                    0 => {
                        if let Some(&lc) = params.first() {
                            p.field_type = FieldType::Isotropic(lc);
                        }
                        String::from(&p.field_type)
                    }
                    1 => {
                        if let [lx, ly, lz] = params[..3] {
                            p.field_type = FieldType::Orthotropic(lx, ly, lz);
                        }
                        String::from(&p.field_type)
                    }
                    _ => "".to_string(),
                };
                let item = StandardListViewItem::from(s.as_str());
                pm.set_row_data(idx as usize, item);
            }
        }
    });

    ui.run()?;

    Ok(())
}
