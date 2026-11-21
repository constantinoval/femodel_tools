use femodel::prelude::*;
use femodel::utils::generations::{FieldType, Perturbation};
use std::io::Write;
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};

fn main() {
    let mut m = FEModel::new();
    m.read("mesh.inp").unwrap();
    let a = Arc::new(AtomicU32::new(0));
    let nnodes = m.nodes.len();
    let h = std::thread::spawn({
        let a = Arc::clone(&a);
        move || {
            Perturbation::new_with_status(&m, FieldType::Isotropic(1.0), a);
        }
    });
    while !h.is_finished() {
        let n = a.load(Ordering::Relaxed);
        print!("\r{:.0}", n as f64 / nnodes as f64 * 100.0);
        std::io::stdout().flush().unwrap();
        // std::thread::sleep(Duration::from_secs(2));
    }
    println!();
    let _hm = h.join().unwrap();
}
