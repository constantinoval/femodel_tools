// Prevent console window in addition to Slint window in Windows release builds when, e.g., starting the app via file manager. Ignored on other platforms.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use femodel::prelude::*;
use rfd::AsyncFileDialog;
use std::error::Error;

slint::include_modules!();

fn main() -> Result<(), Box<dyn Error>> {
    let ui = AppWindow::new()?;
    let ui_handle = ui.as_weak();
    let m = std::rc::Rc::new(std::cell::RefCell::new(FEModel::default()));
    ui.on_choose_file(move || {
        let ui_handle = ui_handle.clone();
        let m = m.clone();
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
                m.borrow_mut().read(&path).unwrap(); // Сохраняем в UI
                println!("{}", m.borrow());
                let m = m.borrow();
                let ncount = m.nodes.len();
                let solids = m.solids.len();
                let shells = m.shells.len();
                let bb = m.bbox();
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
        .unwrap(); // ui.set_model_size("1212".into());
    });
    ui.run()?;

    Ok(())
}
