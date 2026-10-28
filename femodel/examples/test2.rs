// use std::collections::HashMap;

use lsmesh::lsmodel::LsModel;
// use rand::RngExt;
use std::collections::HashMap;

fn main() {
    let mut m = LsModel::default();
    let cell_size = 2.0;
    let radius = 2.0;
    // let mut rng = rand::rng();
    match m.read("shell.k") {
        Ok(()) => {
            let sh = m.nodes_spacial_hash(cell_size);

            // dbg!(&sh[&(2, 0, 0)]);
            // let r: HashMap<usize, f64> = m.nodes.values().map(|n| (n.n, rng.random())).collect();
            let nds = m.node_neibors(&m.nodes[&1], &sh, cell_size, radius);
            dbg!(&nds);
            println!("{} {}", m.nodes.len(), nds.len());
            let mut flag: HashMap<usize, f64> = m.nodes.keys().map(|k| (*k, 0.0)).collect();
            for n in nds {
                flag.insert(n, 1.0);
            }
            m.save_vtk("1.vtk", Some(&[("flag", &flag)]), None).unwrap();
        }
        _ => println!("Error reading file"),
    }
}
