use anyhow::{self, Context};
use femodel::prelude::*;
use femodel::utils::generations::norm_distribution_to_mean_std;
use femodel::utils::generations::{
    FieldProperties, FieldType, GeneratedField, fields_save_to_abaqus,
};
use serde::{Deserialize, Serialize};
use std::env;
use std::fs::File;
use std::io::BufReader;

#[derive(Serialize, Deserialize, Debug)]
struct JsonConfig {
    field_properties: Vec<FieldProperties>,
    coupled: bool,
    mesh_file: String,
    output: String,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    let arg1 = args
        .get(1)
        .context("Необходим параметр - путь к json файлу с конфигурацией генерации")?;
    let f = File::open(arg1)?;
    let reader = BufReader::new(f);
    let config: JsonConfig = serde_json::from_reader(reader)?;
    let mut m = FEModel::default();
    println!("Чтение сетки...");
    m.read(&config.mesh_file).map_err(anyhow::Error::msg)?;
    println!("Генерация случайных полей...");
    let generated_fields = if config.coupled {
        let mut generated_fields: Vec<GeneratedField> = Vec::new();
        let std_normal_field = match config.field_properties[0].field_type {
            FieldType::Isotropic(lc) => m.generate_nodal_correlated_random_field(lc / 2.0),
            FieldType::Orthotropic((lx, ly, lz)) => {
                m.generate_nodal_correlated_random_field_ortho(lx / 2.0, ly / 2.0, lz / 2.0)
            }
        };
        for f in config.field_properties {
            let std = if f.inverse { -f.std } else { f.std };
            let data = norm_distribution_to_mean_std(&std_normal_field, &f.mean, &std);
            generated_fields.push(GeneratedField {
                field_property: f,
                field: data,
            });
        }
        generated_fields
    } else {
        let mut generated_fields: Vec<GeneratedField> = Vec::new();
        for f in config.field_properties {
            println!("Генерация случайного поля {}", f.id);
            let std_normal_field = match f.field_type {
                FieldType::Isotropic(lc) => m.generate_nodal_correlated_random_field(lc / 2.0),
                FieldType::Orthotropic((lx, ly, lz)) => {
                    m.generate_nodal_correlated_random_field_ortho(lx / 2.0, ly / 2.0, lz / 2.0)
                }
            };
            let data = norm_distribution_to_mean_std(&std_normal_field, &f.mean, &f.std);
            generated_fields.push(GeneratedField {
                field_property: f,
                field: data,
            });
        }
        generated_fields
    };
    println!("Сохранение данных...");
    fields_save_to_abaqus(&config.output, &generated_fields)?;
    let names: Vec<String> = generated_fields
        .iter()
        .map(|f| format!("F{}", f.field_property.id))
        .collect();
    let nodal_fields: Vec<NodeScalarField> = names
        .iter()
        .zip(generated_fields.iter())
        .map(|(name, f)| (name.as_str(), &f.field))
        .collect();
    let file_path = std::path::Path::new(&config.output);
    let vtk_path = file_path.with_extension("vtk");
    m.save_vtk(vtk_path.to_str().unwrap(), Some(&nodal_fields), None)?;
    Ok(())
}
