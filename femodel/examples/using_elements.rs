use lsmesh::lsmodel::lselement::LsElement;

fn main() {
    let mut e: LsElement = "1, 1, 1, 2, 3, 4, 5, 5, 6, 6".into();
    e.normalize_nodes();
    println!("{:?}", e.nodes);
    println!("{}", e.nodes_count());
    println!("{:?}", &e);
    println!("{}", String::try_from(&e).unwrap());
    let mut e: LsElement = "1, 1, 1, 2, 3".into();
    e.etype = lsmesh::lsmodel::lselement::ElementType::SHELL;
    e.normalize_nodes();
    println!("{}", String::try_from(&e).unwrap());
    println!("{}", String::try_from(&e).unwrap());
}
