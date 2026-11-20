// Prevent console window in addition to Slint window in Windows release builds when, e.g., starting the app via file manager. Ignored on other platforms.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use femodel::prelude::*;
use femodel::utils::generations::{FieldType, Perturbation};
use rfd::AsyncFileDialog;
use slint::{Model, ModelRc, StandardListViewItem, VecModel};
use std::error::Error;
use std::rc::Rc;

struct FieldDraft {
    id: usize,
    name: String,
    mean: f32,
    std: f32,
    inverse: bool,
}

struct PerturbationDraft {
    field_type: FieldType,
    fields: Vec<FieldDraft>,
}

struct AppState {
    model: FEModel, // None пока файл не открыт
    perturbations: Vec<PerturbationDraft>,
}

slint::include_modules!();

fn main() -> Result<(), Box<dyn Error>> {
    let ui = AppWindow::new()?;
    let ui_handle = ui.as_weak();
    let app_state = Rc::new(std::cell::RefCell::new(AppState {
        model: FEModel::default(),
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
                    app_state.borrow_mut().model.read(&path).unwrap(); // Сохраняем в UI
                    println!("{}", app_state.borrow().model);
                    let ncount = app_state.borrow().model.nodes.len();
                    let solids = app_state.borrow().model.solids.len();
                    let shells = app_state.borrow().model.shells.len();
                    let bb = app_state.borrow().model.bbox();
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
        move || {
            let app_state = app_state.clone();
            slint::spawn_local(async move {
                let app_state = app_state.borrow();
                let pert =
                    Perturbation::new(&app_state.model, app_state.perturbations[0].field_type);
            })
            .unwrap();
        }
    });

    ui.run()?;

    Ok(())
}
