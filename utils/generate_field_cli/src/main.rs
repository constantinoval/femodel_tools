use femodel::prelude::*;
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

fn ask_positive_float(message: &str, dflt: f64) -> f64 {
    Text::new(message)
        .with_validator(|input: &str| {
            // поддержка запятой как разделителя
            let normalized = input.trim().replace(',', ".");
            match normalized.parse::<f64>() {
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
        .parse::<f64>()
        .unwrap()
}

fn ask_nonnegative_float(message: &str, dflt: f64) -> f64 {
    Text::new(message)
        .with_validator(|input: &str| {
            // поддержка запятой как разделителя
            let normalized = input.trim().replace(',', ".");
            match normalized.parse::<f64>() {
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
        .parse::<f64>()
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

#[derive(Debug)]
enum FieldType {
    Isotropic(f64),
    Orthotropic((f64, f64, f64)),
}

#[derive(Debug)]
struct FieldProperties {
    mean: f64,
    std: f64,
    id: u16,
    name: String,
    field_type: FieldType,
}

impl FieldProperties {
    fn new() -> Self {
        Self {
            mean: 0.0,
            std: 1.0,
            id: 1,
            name: "".to_string(),
            field_type: FieldType::Isotropic(1.0),
        }
    }
}

fn main() -> anyhow::Result<()> {
    let files = find_files(&["inp", "k"]);
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

    let mut fields: Vec<FieldProperties> = Vec::new();

    loop {
        let mut fp = FieldProperties::new();
        fp.mean = ask_nonnegative_float("Среднее значение генерируемого поля:", 0.0);
        fp.std = ask_positive_float("Стандартное отклонение генерируемого поля:", 1.0);
        fp.id = ask_positive_int("Номер поля для записи в файл:", 1);
        fp.name = Text::new("Имя или комментарий к полю:")
            .with_help_message("Например: напряжение течения | ef")
            .with_default("")
            .prompt()
            .expect("Ошибка ввода");
        let field_type = Select::new("Тип неоднородности:", vec!["изотропная", "ортотропная"])
            .prompt()
            .expect("Ошибка выбора");
        match field_type {
            "изотропная" => {
                let lc = ask_positive_float("Масштаб корреляции Lc:", 1.0);
                fp.field_type = FieldType::Isotropic(lc);
            }
            "ортотропная" => {
                let lx = ask_positive_float("Масштаб корреляции по оси 0X Lx:", 1.0);
                let ly = ask_positive_float("Масштаб корреляции по оси 0Y Ly:", 1.0);
                let lz = ask_positive_float("Масштаб корреляции по оси 0Z Lz:", 1.0);
                fp.field_type = FieldType::Orthotropic((lx, ly, lz));
            }
            _ => {}
        }
        fields.push(fp);
        if !Confirm::new("Добавить еще поля?")
            .with_default(false)
            .prompt()
            .expect("Ошибка ввода")
        {
            break;
        }
    }
    let f = m.generate_nodal_correlated_random_field(1.0);
    dbg!(fields);
    let coupled_fields =
        Confirm::new("Генерировать поля с одинаковым распределением в пространстве?")
            .with_default(false)
            .prompt()
            .expect("Ошибка ввода");
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
    // let Lc = ask_non_negative_text("Введите размер корреляции:");
    Ok(())
}
