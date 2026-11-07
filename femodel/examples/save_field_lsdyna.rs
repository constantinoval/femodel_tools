use anyhow;
use femodel::prelude::*;
use femodel::utils::generations::{FieldProperties, FieldType, GeneratedField};

fn main() -> anyhow::Result<()> {
    let mut m = FEModel::new();
    m.read("mesh1.k").expect("Can't read model");
    let fp = FieldProperties {
        mean: 0.02,
        std: 0.005,
        id: 6,
        name: "D".to_string(),
        field_type: FieldType::Isotropic(10.0),
        inverse: false,
    };
    let field = m.generate_field_using_property_object(&fp);
    let gf = vec![GeneratedField {
        field_property: fp,
        field,
    }];
    femodel::utils::generations::save_fields_to_lsdyna("damage.k", &gf, true)?;
    Ok(())
}
