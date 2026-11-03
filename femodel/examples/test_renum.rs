use femodel::prelude::*;
use std::time::Instant;

fn main() {
    let start = Instant::now();
    let mut m = FEModel::default();
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
            let _ = m.save("1.k");
        }
        Err(msg) => println!("Error: {msg}"),
    }
}
