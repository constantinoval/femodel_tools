use lsmesh::lsmodel::LsModel;
use std::time::Instant;

fn main() {
    let mut m = LsModel::default();
    let start = Instant::now();
    match m.read("last_state_3D.k") {
        Ok(()) => {
            let duration = start.elapsed();
            println!(
                "num of nodes {}, num of solids {}, num of shells {}, duration: {:?}",
                m.nodes.len(),
                m.solids.len(),
                m.shells.len(),
                duration
            );
        }
        Err(msg) => println!("Error: {msg}"),
    }
    println!("Searching for unconnected regions...");
    let t0 = Instant::now();
    let c = m.find_unconnected_regions(&[]);
    println!(
        "{} unconnected regions, search duration: {:?}",
        c.len(),
        t0.elapsed()
    );
}
