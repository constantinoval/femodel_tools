use lsmesh::lsmodel::LsModel;
use lsmesh::lsmodel::lsnode::LsNode;
use lsmesh::lsmodel::utils::renumber_map;
use std::collections::{BTreeMap, HashSet};
use std::time::Instant;

fn main() {
    let start = Instant::now();
    let mut m = LsModel::default();
    match m.read("2_shells_2_solids.k") {
        Ok(()) => {
            let duration = start.elapsed();
            println!(
                "num of nodes {}, num of solids {}, num of shells {}, duration: {:?}",
                m.nodes.len(),
                m.solids.len(),
                m.shells.len(),
                duration
            );
            println!("{:?}, {:?}", m.solids, m.shells);
            m.renumber_elements(100);
            println!("{:?}, {:?}", m.solids, m.shells);
            m.renumber_nodes(1000);
            println!("{:?}, {:?}", m.solids, m.shells);
            m.renumber_parts(10);
            println!("{:?}, {:?}", m.solids, m.shells);
            m.save("1.k");
        }
        Err(msg) => println!("Error: {msg}"),
    }
}
