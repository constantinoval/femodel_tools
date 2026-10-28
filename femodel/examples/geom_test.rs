use lsmesh::lsmodel::geom_utils::Vector3D;
use lsmesh::vec3d;

fn main() {
    let mut p = Vector3D {
        x: 2.0,
        y: 0.0,
        z: 0.0,
    };
    println!("Point: {:?}", p);
    println!("Point norm: {}", p.norm());
    println!("After normalization");
    p.normalize();
    assert_eq!(p.norm(), 1.0f64);
    println!("Point: {:?}", p);
    println!("Point norm: {}", p.norm());
    p.x = 4.0;
    let p2 = p.normalized();
    println!("Point 2: {:?}", p2);
    println!("Point 2 norm: {}", p2.norm());
    let p = vec3d![13.0, 0.0, 0.0];
    println!("{:?}", p);
    let p2 = vec3d![1.0, 3.0, 0.0];
    let s = p.dot(&p2);
    println!("Scalar product of the vectors is {}", s);
    let mut p = vec3d![1.0, 0.0, 0.0];
    p.scale(10);
    println!("{:?}", p);
}
