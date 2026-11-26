use anyhow::{self, Context};
use femodel::prelude::*;
use femodel::utils::generations::{
    FieldType, Perturbation, RandomField, save_fields_to_abaqus, save_fields_to_lsdyna,
};
use serde::{Deserialize, Serialize};
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

#[derive(Serialize, Deserialize, Debug)]
struct FieldPropertiesJson {
    id: usize,
    name: String,
    mean: f32,
    std: f32,
    invert: bool,
    trunc: Option<(f32, f32)>,
}

#[derive(Serialize, Deserialize, Debug)]
struct RandomFieldsJson {
    field_type: FieldType,
    properties: Vec<FieldPropertiesJson>,
}

#[derive(Serialize, Deserialize, Debug)]
struct JsonConfig {
    mesh_file: String,
    output: String,
    save_vtk: bool,
    fields: Vec<RandomFieldsJson>,
}

fn save_json_example(file_path: &str) -> anyhow::Result<()> {
    let fp1 = FieldPropertiesJson {
        id: 1,
        name: "field_1".into(),
        mean: 100.0,
        std: 10.0,
        invert: false,
        trunc: None,
    };
    let fp2 = FieldPropertiesJson {
        id: 2,
        name: "field_2".into(),
        mean: 50.0,
        std: 2.0,
        invert: true,
        trunc: Some((48.0, 52.0)),
    };
    let rf1 = RandomFieldsJson {
        field_type: FieldType::Isotropic(1.0),
        properties: vec![fp1, fp2],
    };
    let fp3 = FieldPropertiesJson {
        id: 3,
        name: "field_3".into(),
        mean: 1.0,
        std: 0.1,
        invert: false,
        trunc: None,
    };
    let rf2 = RandomFieldsJson {
        field_type: FieldType::Orthotropic(1.0, 2.0, 3.0),
        properties: vec![fp3],
    };
    let jc = JsonConfig {
        mesh_file: "mesh.k".into(),
        output: "fieds.k".into(),
        save_vtk: true,
        fields: vec![rf1, rf2],
    };

    let s = serde_json::to_string_pretty(&jc)?;
    std::fs::write(file_path, s)?;
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    println!("Генерация случайных коррелированных полей по данным из конфигурационного файла.");
    println!("Вример вызова: generate_field_json config.json");
    println!(
        "Если указанный конфигурационный файл не существует, будет записан пример конфигурационного файла."
    );

    let arg1 = args
        .get(1)
        .context("Необходим параметр - путь к json файлу с конфигурацией генерации")?;
    let pb = PathBuf::from(&arg1);
    if !pb.exists() {
        println!("Файл {arg1} не найден. Запись примера конфигурационного файла в {arg1}");
        save_json_example(&arg1)?;
        return Ok(());
    }
    let f = File::open(arg1)?;
    let reader = BufReader::new(f);
    let config: JsonConfig = serde_json::from_reader(reader)?;
    let mut m = FEModel::default();
    println!("Чтение сетки...");
    m.read(&config.mesh_file).map_err(anyhow::Error::msg)?;
    let mut random_fields: Vec<RandomField> = Vec::new();
    println!("Генерация случайных полей...");
    for fields in config.fields {
        let pert = Perturbation::new(&m, fields.field_type);
        for props in fields.properties {
            let rf = RandomField::new()
                .with_id(props.id)
                .with_name(props.name)
                .generate(&pert, props.mean, props.std, props.invert, props.trunc);
            random_fields.push(rf);
        }
    }
    println!("Сохранение данных...");
    let file_path = std::path::Path::new(&config.output);
    match file_path.extension().and_then(|e| e.to_str()) {
        Some("inp") => save_fields_to_abaqus(
            &config.output,
            &random_fields
                .iter()
                .map(|rf| rf)
                .collect::<Vec<&RandomField>>(),
        )?,
        Some("k") => save_fields_to_lsdyna(
            &config.output,
            &random_fields
                .iter()
                .map(|rf| rf)
                .collect::<Vec<&RandomField>>(),
        )?,
        Some(other) => println!("File type with extension {other} not supported"),
        None => println!("There is no extension in output file"),
    }
    if config.save_vtk {
        let names: Vec<String> = random_fields.iter().map(|f| f.name.clone()).collect();
        let nodal_fields: Vec<NodeScalarField> = names
            .iter()
            .zip(random_fields.iter())
            .map(|(name, f)| (name.as_str(), &f.data))
            .collect();
        let vtk_path = file_path.with_extension("vtk");
        m.save_vtk(vtk_path.to_str().unwrap(), Some(&nodal_fields), None)?;
    }
    Ok(())
}
