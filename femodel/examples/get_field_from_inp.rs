use femodel::readers::abaqus::read_nodal_fields_from_inp;
use std::time::Instant;

fn main() {
    let t0 = Instant::now();
    let f = read_nodal_fields_from_inp("1.inp").unwrap();
    dbg!(t0.elapsed());
    dbg!(f.keys());
    println!("Len of data {}", f[&1].len());
}
