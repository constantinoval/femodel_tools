use lsmesh::lsmodel::LsModel;
use rayon::prelude::*;
use std::time::Instant;

fn main() {
    let mut m = LsModel::default();
    let start = Instant::now();
    match m.read("mesh.k") {
        Ok(()) => {
            let duration = start.elapsed();
            println!(
                "num of nodes {}, num of solids {}, num of shells {}, duration: {:?}",
                m.nodes.len(),
                m.solids.len(),
                m.shells.len(),
                duration
            );
            println!("{:?}, {}", m.solids[&1].nodes, m.solids[&1].nodes_count());
        }
        Err(msg) => println!("Error: {msg}"),
    }
    println!("Calculating volumes");
    let start = Instant::now();
    m.calculate_element_volumes();
    let duration = start.elapsed();
    println!("Done in {:?}", duration);
    println!("Volume of element 1 is {}", m.solids[&1].volume);
    let v: f64 = m.solids.par_iter().map(|(_, e)| e.volume).sum();
    println!("Total volume is {}", v);
    let start = Instant::now();
    println!(
        "There are {} elements in part 1 (check duration: {:?})",
        m.elements_by_part(&1).len(),
        start.elapsed()
    );
    let start = Instant::now();
    println!(
        "There are {} elements in part 2 (check duration: {:?})",
        m.elements_by_part(&2).len(),
        start.elapsed()
    );
}
