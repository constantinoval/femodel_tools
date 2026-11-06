use femodel::prelude::*;

fn main() {
    let mut m1 = FEModel::new();
    m1.read("mesh2.k").unwrap();
    let mut m2 = FEModel::new();
    m2.read("mesh1.k").unwrap();
    let f1 = m1.generate_nodal_correlated_random_field(10.0);
    // let f2 = m1.remap_to_model_by_nearest(10.0, &m2, &f1);
    let f2 = m1.remap_to_model_by_idw(10.0, &m2, &f1);
    m1.save_vtk("m1.vtk", Some(&[("f1", &f1)]), None).unwrap();
    m2.save_vtk("m2.vtk", Some(&[("f1", &f2)]), None).unwrap();
}
