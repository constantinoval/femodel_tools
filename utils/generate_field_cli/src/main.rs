use femodel::objects::femodel::NodeScalarField;
use femodel::prelude::*;
use femodel::utils::generations::{
    FieldType, Perturbation, RandomField, save_fields_to_abaqus, save_fields_to_lsdyna,
};
use inquire::{Confirm, Select, Text, validator::Validation};
use std::fs;

fn find_files(exts: &[&str]) -> Vec<String> {
    fs::read_dir(".")
        .expect("Не удалось прочитать директорию")
        .filter_map(|entry| entry.ok())
        .filter(|e| e.file_type().map(|t| t.is_file()).unwrap_or(false))
        .filter_map(|e| {
            let path = e.path();
            let ext = path.extension().and_then(|s| s.to_str())?;
            if exts.iter().any(|x| x.eq_ignore_ascii_case(ext)) {
                path.file_name()
                    .and_then(|n| n.to_str())
                    .map(|s| s.to_string())
            } else {
                None
            }
        })
        .collect()
}

fn ask_positive_float(message: &str, dflt: f32) -> f32 {
    Text::new(message)
        .with_validator(|input: &str| {
            // поддержка запятой как разделителя
            let normalized = input.trim().replace(',', ".");
            match normalized.parse::<f32>() {
                Ok(v) if v <= 0.0 => Ok(Validation::Invalid(
                    "Число должно быть положительным".into(),
                )),
                Ok(_) => Ok(Validation::Valid),
                Err(_) => Ok(Validation::Invalid("Введите корректное число".into())),
            }
        })
        .with_default(format!("{}", dflt).as_str())
        .prompt()
        .expect("Ошибка ввода")
        .trim()
        .replace(',', ".")
        .parse::<f32>()
        .unwrap()
}

fn ask_nonnegative_float(message: &str, dflt: f32) -> f32 {
    Text::new(message)
        .with_validator(|input: &str| {
            // поддержка запятой как разделителя
            let normalized = input.trim().replace(',', ".");
            match normalized.parse::<f32>() {
                Ok(v) if v < 0.0 => Ok(Validation::Invalid(
                    "Число должно быть неотрцательным".into(),
                )),
                Ok(_) => Ok(Validation::Valid),
                Err(_) => Ok(Validation::Invalid("Введите корректное число".into())),
            }
        })
        .with_default(format!("{}", dflt).as_str())
        .prompt()
        .expect("Ошибка ввода")
        .trim()
        .replace(',', ".")
        .parse::<f32>()
        .unwrap()
}

fn ask_positive_int(message: &str, dflt: u16) -> u16 {
    Text::new(message)
        .with_default(format!("{}", dflt).as_str())
        .prompt()
        .expect("Ошибка ввода")
        .trim()
        .parse::<u16>()
        .unwrap()
}

fn main() -> anyhow::Result<()> {
    let files = find_files(&["inp", "k"]);
    let mut field_num: u16 = 1;
    if files.is_empty() {
        anyhow::bail!(
            "В текущей директории не найдено файлов с сеточной моделью Abaqus (inp) или ls-dyna (k)"
        );
    }
    let file = Select::new("Выберите файл:", files)
        .prompt()
        .expect("Ошибка выбора");

    let mut m = FEModel::default();
    println!("Чтение файла {}", file);
    m.read(&file).map_err(anyhow::Error::msg)?;
    let bb = m.bbox();
    println!(
        "Габариты модели: x: {} -> {}, y: {}->{}, z: {}->{}",
        bb.min_x, bb.max_x, bb.min_y, bb.max_y, bb.min_z, bb.max_z
    );

    let mut fields: Vec<RandomField> = Vec::new();

    let num_perturbations = ask_positive_int("Число распределений? ", 1);
    for n in 0..num_perturbations {
        println!("Наполнение распределения {}...", n + 1);
        let field_type = Select::new("Тип неоднородности:", vec!["изотропная", "ортотропная"])
            .prompt()
            .expect("Ошибка выбора");
        let field_type: FieldType = match field_type {
            "изотропная" => {
                let lc = ask_positive_float("Масштаб корреляции Lc:", 1.0);
                Some(FieldType::Isotropic(lc))
            }
            "ортотропная" => {
                let lx = ask_positive_float("Масштаб корреляции по оси 0X Lx:", 1.0);
                let ly = ask_positive_float("Масштаб корреляции по оси 0Y Ly:", 1.0);
                let lz = ask_positive_float("Масштаб корреляции по оси 0Z Lz:", 1.0);
                Some(FieldType::Orthotropic(lx, ly, lz))
            }
            _ => None,
        }
        .unwrap();
        let pert = Perturbation::new(&m, field_type);
        println!("Добавляем поля к распределению {}...", n + 1);
        loop {
            let id = ask_positive_int("Номер поля для записи в файл:", field_num);
            field_num = id + 1;
            let name = Text::new("Имя или комментарий к полю:")
                .with_help_message("Например: напряжение течения | ef")
                .with_default("")
                .prompt()
                .expect("Ошибка ввода");
            let mean = ask_nonnegative_float("Среднее значение генерируемого поля:", 0.0);
            let std = ask_positive_float("Стандартное отклонение генерируемого поля:", 1.0);
            let inverse = Confirm::new("Инвертировать поле?")
                .with_default(false)
                .prompt()
                .expect("Ошибка ввода");
            let field = RandomField::new()
                .with_id(id as usize)
                .with_name(name)
                .generate(&pert, mean, std, inverse);
            fields.push(field);
            if !Confirm::new("Добавить еще поля?")
                .with_default(false)
                .prompt()
                .expect("Ошибка ввода")
            {
                break;
            }
        }
    }
    let path = std::path::Path::new(&file);
    let name = path
        .file_stem()
        .and_then(|f| f.to_str())
        .unwrap_or("")
        .to_string();
    let ext = path
        .extension()
        .and_then(|f| f.to_str())
        .unwrap()
        .to_string();
    let file_name = Text::new("Имя файла для сохранения:")
        .with_default(format!("{}_{}.{}", name, "fields", ext).as_str())
        .prompt()
        .expect("Ошибка ввода");
    let file_path = std::path::Path::new(&file_name);
    match file_path.extension().and_then(|e| e.to_str()) {
        Some("inp") => save_fields_to_abaqus(
            &file_name,
            &fields.iter().map(|rf| rf).collect::<Vec<&RandomField>>(),
        )?,
        Some("k") => save_fields_to_lsdyna(
            &file_name,
            &fields.iter().map(|rf| rf).collect::<Vec<&RandomField>>(),
        )?,
        Some(other) => println!("File type with extension {other} not supported"),
        None => println!("There is no extension in output file"),
    }
    if Confirm::new("Записать vtk? ")
        .with_default(false)
        .prompt()
        .expect("Ошибка ввода")
    {
        let names: Vec<String> = fields.iter().map(|f| f.name.clone()).collect();
        let nodal_fields: Vec<NodeScalarField> = names
            .iter()
            .zip(fields.iter())
            .map(|(name, f)| (name.as_str(), &f.data))
            .collect();
        let vtk_path = file_path.with_extension("vtk");
        m.save_vtk(vtk_path.to_str().unwrap(), Some(&nodal_fields), None)?;
    }
    Ok(())
}
