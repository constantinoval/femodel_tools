// Prevent console window in addition to Slint window in Windows release builds when, e.g., starting the app via file manager. Ignored on other platforms.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use rfd::FileDialog;
use std::error::Error;

slint::include_modules!();

fn main() -> Result<(), Box<dyn Error>> {
    let ui = AppWindow::new()?;
    let ui_handle = ui.as_weak();

    ui.on_choose_file(move || {
        let ui = ui_handle.unwrap();
        let pick = FileDialog::new()
            .add_filter("supported meshes", &["k", "inp"])
            .add_filter("ls-dyna mesh", &["k"])
            .add_filter("abaqus mesh", &["inp"])
            .pick_file();
        // ui.set_model_size("1212".into());
        let mesh_path = if let Some(path) = pick {
            path.to_str().unwrap().to_owned().into()
        } else {
            "".into()
        };
        mesh_path
    });
    ui.run()?;

    Ok(())
}
