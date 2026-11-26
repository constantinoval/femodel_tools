use femodel::prelude::*;
use femodel::utils::generations;

fn main() -> Result<(), String> {
    let mut m = FEModel::default();
    m.read("mesh.inp")?;
    let q = m.generate_nodal_correlated_random_field(1.5);
    // let q = generations::norm_to_damage_lognormal(&q, 0.01, 0.2);
    // let q = generations::norm_to_uniform_rank(&q);
    // let q = generations::norm_trunc(&q, 0.0, 1.0);
    let q = generations::transform_to_truncated_normal(&q, 0.0, 0.1, 0.1);
    m.save_vtk("1.vtk", Some(&[("d", &q)]), None).unwrap();
    Ok(())
}
