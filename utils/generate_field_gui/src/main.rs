// Prevent console window in addition to Slint window in Windows release builds when, e.g., starting the app via file manager. Ignored on other platforms.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use femodel::prelude::*;
use femodel::utils::generations::FieldType;
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
    let app = app_state.clone();
    let ui_handle1 = ui_handle.clone();
    ui.on_choose_file(move || {
        let ui_handle = ui_handle1.clone();
        let app_state = app.clone();
        // let mut m = &app_state.model;
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
                let mut s =
                    format!("Узлов: {ncount}, объемных элементов: {solids}, оболочек: {shells}\n");
                s += &format!(
                    "Габариты модели: x {}->{}, y: {}->{}, z: {}->{}",
                    bb.min_x, bb.max_x, bb.min_y, bb.max_y, bb.min_z, bb.max_z
                );
                let _ = ui_handle.upgrade_in_event_loop(move |ui| {
                    ui.set_file_path(path.into());
                    ui.set_model_info(s.into());
                });
            }
        })
        .unwrap();
    });

    let app = app_state.clone();
    let pm = perturbation_model.clone();
    ui.on_add_perturbation(move |field_type: i32, params: ModelRc<f32>| {
        let app_state = app.clone();
        let params: Vec<f32> = params.iter().collect();
        let ui = ui_handle.unwrap();
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
                    String::from(&FieldType::Orthotropic(params[0], params[1], params[2])).as_str(),
                ));
            }
            _ => {}
        }
        let idx = (app_state.borrow().perturbations.len() - 1) as i32;
        ui.set_current_perturbation_index(idx);
    });

    let app = app_state.clone();
    let pm = perturbation_model.clone();
    ui.on_delete_current_perturbation(move |idx: i32| {
        if idx == -1 {
            return;
        }
        let app_state = app.clone();
        // let perturbation_model = pm;
        app_state.borrow_mut().perturbations.remove(idx as usize);
        pm.remove(idx as usize);
    });

    ui.run()?;

    Ok(())
}
