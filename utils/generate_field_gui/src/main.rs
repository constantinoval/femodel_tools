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
use std::path::PathBuf;
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
        move |idx, id, name, mean, std, inverse| {
            app_state.borrow_mut().perturbations[idx as usize]
                .fields
                .push(FieldDraft {
                    id: id as usize,
                    name: String::from(&name),
                    mean,
                    std,
                    inverse,
                });
            fm.push(StandardListViewItem::from(
                String::from(format!(
                    "id({}), name({}), mean({}), std({}), inverse({})",
                    id, name, mean, std, inverse
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
            let save_vtk = ui.get_is_save_vtk();
            let st = app_state.borrow();
            let m = st.model.read().unwrap();
            if m.nodes.len() == 0 {
                return;
            }
            if st.perturbations.len() == 0 {
                return;
            }
            let app_state = app_state.clone();
            let ui_handle = ui_handle.clone();

            // Блокируем UI сразу — до диалога.
            ui_handle
                .upgrade_in_event_loop(|ui| ui.set_generating(true))
                .ok();

            slint::spawn_local(async move {
                let reset = {
                    let ui_handle = ui_handle.clone();
                    move || {
                        ui_handle
                            .upgrade_in_event_loop(|ui| ui.set_generating(false))
                            .ok();
                    }
                };

                let pick = AsyncFileDialog::new()
                    .add_filter("ls-dyna", &["k"])
                    .add_filter("abaqus", &["inp"])
                    .set_file_name("fields")
                    .save_file()
                    .await;
                let Some(file) = pick else {
                    reset();
                    return;
                };
                let out_path = file.path().to_path_buf();

                let (model_arc, jobs) = {
                    let st = app_state.borrow();
                    (Arc::clone(&st.model), st.perturbations.clone())
                };
                if jobs.iter().all(|p| p.fields.is_empty()) {
                    eprintln!("Нет полей для генерации");
                    reset();
                    return;
                }

                let ui_handle = ui_handle.clone();
                std::thread::spawn(move || {
                    let model = model_arc.read().unwrap();
                    let results: Vec<RandomField> = jobs
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
                    let ext = out_path
                        .extension()
                        .and_then(|e| e.to_str())
                        .unwrap_or("")
                        .to_ascii_lowercase();
                    let _ = match ext.as_str() {
                        "k" => save_fields_to_lsdyna(out_path.to_str().unwrap(), &refs),
                        _ => save_fields_to_abaqus(out_path.to_str().unwrap(), &refs),
                    };
                    if save_vtk {
                        let names: Vec<String> = results.iter().map(|f| f.name.clone()).collect();
                        let nodal_fields: Vec<NodeScalarField> = names
                            .iter()
                            .zip(results.iter())
                            .map(|(name, f)| (name.as_str(), &f.data))
                            .collect();
                        let vtk_path = out_path.with_extension("vtk");
                        let _ = model.save_vtk(
                            vtk_path.to_str().unwrap_or_default(),
                            Some(&nodal_fields),
                            None,
                        );
                    }
                    drop(model);
                    let _ = ui_handle.upgrade_in_event_loop(move |ui| {
                        ui.set_generating(false);
                    });
                });
            })
            .ok(); // не unwrap: без паники, если очередь недоступна
        }
    });

    ui.run()?;

    Ok(())
}
