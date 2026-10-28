use std::collections::HashMap;

use lsmesh::lsmodel::LsModel;
use rand::RngExt;
// use std::collections::HashMap;

fn main() {
    let mut m = LsModel::default();
    let mut rng = rand::rng();
    match m.read("27_solids.k") {
        Ok(()) => {
            // let sh = m.spacial_hash(0.1);
            // dbg!(sh);
            // println!("{:?}", sh.get(&(0, 5, 0)));
            // let nn = sh.values().fold(0, |acc, v| acc + v.len());
            // println!(
            //     "Nodes count: {}, nodes count from sh: {}",
            //     m.nodes.len(),
            //     nn
            // );
            let r: HashMap<usize, f64> = m.nodes.values().map(|n| (n.n, rng.random())).collect();
            let x: HashMap<usize, f64> = m.nodes.values().map(|n| (n.n, n.x)).collect();
            let y: HashMap<usize, f64> = m.nodes.values().map(|n| (n.n, n.y)).collect();
            let z: HashMap<usize, f64> = m.nodes.values().map(|n| (n.n, n.z)).collect();
            // let n: HashMap<usize, f64> = m.nodes.values().map(|n| (n.n, n.n as f64)).collect();
            // let i: HashMap<usize, f64> = m
            //     .nodes
            //     .values()
            //     .map(|n| (n.n, n.spacial_hash(0.1).0 as f64))
            //     .collect();
            // let j: HashMap<usize, f64> = m
            //     .nodes
            //     .values()
            //     .map(|n| (n.n, n.spacial_hash(0.1).1 as f64))
            //     .collect();
            // let k: HashMap<usize, f64> = m
            //     .nodes
            //     .values()
            //     .map(|n| (n.n, n.spacial_hash(0.1).2 as f64))
            //     .collect();
            // let e_n: HashMap<usize, f64> = m.solids.values().map(|e| (e.n, e.n as f64)).collect();
            // m.save_vtk(
            //     "1.vtk",
            //     Some(&[
            //         ("x", &x),
            //         ("y", &y),
            //         ("z", &z),
            //         ("n", &n),
            //         ("i", &i),
            //         ("j", &j),
            //         ("k", &k),
            //     ]),
            //     Some(&[("e_n", &e_n)]),
            // )
            // .unwrap();
            m.save_vtk("1.vtk", Some(&[("x", &x), ("y", &y), ("z", &z)]), None)
                .unwrap();
        }
        _ => println!("Error reading file"),
    }
}
