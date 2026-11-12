use femodel::prelude::*;
use femodel::utils::generations::{
    FieldType, Perturbation, RandomField, save_fields_to_abaqus, save_fields_to_lsdyna,
};

fn main() {
    let mut m = FEModel::new();
    m.read("mesh.k").unwrap();
    let p = Perturbation::new(&m, FieldType::Isotropic(10.0));
    for n in p.node_ids() {
        println!("id {n}");
    }
    let v = p.node_normalized_value(1, 100.0, 10.0).unwrap();
    println!("v: {v}");
    let v: Vec<f32> = p
        .node_ids()
        .filter_map(|&n| p.node_normalized_value(n, 100.0, 10.0).ok())
        .collect();
    dbg!(v);
    let v = p.generate_for_norm_std(100.0, 10.0);
    dbg!(v);
    let rf1 = RandomField::new()
        .with_id(1)
        .with_name("s0".to_string())
        .generate(&p, 100.0, 20.0, false);
    let rf2 = RandomField::new()
        .with_id(2)
        .with_name("Et".to_string())
        .generate(&p, 1000.0, 100.0, false);
    let rf3 = RandomField::new()
        .with_id(3)
        .with_name("ef".to_string())
        .generate(&p, 0.5, 0.01, true);
    save_fields_to_lsdyna("1.k", &vec![&rf1, &rf2, &rf3]).unwrap();
    save_fields_to_abaqus("1.inp", &vec![&rf1, &rf2, &rf3]).unwrap();
}
