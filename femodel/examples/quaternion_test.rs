use lsmesh::lsmodel::geom_utils::Quarternion;
use lsmesh::lsmodel::geom_utils::Vector3D;
use lsmesh::vec3d;

fn main() {
    let q = Quarternion::new(180.0, vec3d![0.0, 10.0, 0.0]);
    println!("{:?}", q);
    let mut p = vec3d![15.0, 0.0, 0.0];
    println!("{:?}", p);
    q.rotate_point(&mut p);
    println!("{:?}", p);
    let p2 = q.rotated_point(&p);
    println!("{:?}", p2);
}
