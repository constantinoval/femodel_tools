use femodel::prelude::*;

fn main() {
    let mut n = FENode {
        x: 0.0,
        y: 0.0,
        z: 0.0,
        n: 1,
    };
    println!("{}", n);
    n.translate(&[1.0, 0.0, 0.0]);
    println!("{}", n);
    n.translate(&[0.0, 1.0, 0.0]);
    println!("{}", n);
    n.translate(&[]);
    println!("{}", n);
    println!("Node x is {}", n.x);
    let mut crds = n.coordinates();
    println!("{:?}", crds);
    crds[1] = &2.0;
    println!("{:?}", crds);
    println!("{}", n);
    let crds = n.coordinates_mut();
    println!("{:?}", crds);
    *crds[1] = 2.0;
    println!("{:?}", crds);
    println!("{}", n);
    n.x = 10.0;
    println!("{}", n);
    // let n: u32 = parse("       100", 0, 10);
    // println!("{}", n);
    let n: FENode = "     100            10.1            20.2            30.3".into();
    println!("{}", n);
    let n: FENode = "1, 2.1, 3.2\n".into();
    dbg!(&n);
    // let s: String = (&n).into();
    println!("{}", <String>::from(&n));
    let n = FENode::from("       1               2               0              -5");
    println!("{:?}", n.spacial_hash(1.0));
}
