use femodel::utils::generations::{FieldProperties, FieldType};

fn main() {
    println!("Hello, world!");
    let mut f = FieldProperties::new();
    f.field_type = FieldType::Orthotropic((1.0, 2.0, 3.0));
    let s = serde_json::to_string(&f).unwrap();
    println!("{}", s);
}
