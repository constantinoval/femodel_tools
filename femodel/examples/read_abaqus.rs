use lsmesh::lsmodel::read_flat_abaqus_input;

fn main() {
    let (nodes, solids, shells) = read_flat_abaqus_input("mesh.inp").unwrap();
    println!("{}", nodes.len());
    dbg!(&solids[&1]);
}
