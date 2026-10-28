use anyhow::Result;
use lsmesh::lsmodel::LsModel;

fn main() -> Result<()> {
    let mut m = LsModel::default();
    println!("Reading mesh");
    m.read("mesh.inp").map_err(anyhow::Error::msg)?;
    println!("Generating random field");
    let z = m.generate_nodal_correlated_random_field(1.0);
    println!("Converting to cell values...");
    let v = z.iter().map(|(_, &v)| v).collect::<Vec<f32>>();
    let (mean, std) = lsmesh::lsmodel::utils::mean_and_std(&v);
    println!("Nodal mean: {}, std: {}", mean, std);
    let cell_z = m.nodal_field_to_cell_field(&z);
    // let cell_z = lsmesh::lsmodel::utils::standardize_field(&cell_z);
    let v = cell_z.iter().map(|(_, &v)| v).collect::<Vec<f32>>();
    let (mean, std) = lsmesh::lsmodel::utils::mean_and_std(&v);
    println!("Cell mean: {}, std: {}", mean, std);

    // calculating mean and std
    println!("Calculating ef...");
    let ef = lsmesh::lsmodel::utils::norm_distribution_to_mean_std(&cell_z, 0.6, 0.02);
    let v = ef.iter().map(|(_, &v)| v).collect::<Vec<f32>>();
    let (mean, std) = lsmesh::lsmodel::utils::mean_and_std(&v);
    println!("ef mean: {}, std: {}", mean, std);
    println!("Saving vtk");
    m.save_vtk(
        "1.vtk",
        Some(&[("Z", &z)]),
        Some(&[("ze", &cell_z), ("ef", &ef)]),
    )?;
    println!("All done...");
    Ok(())
}
