use crate::objects::feelement::ElementType;
use crate::prelude::*;
use crate::readers::common_types::*;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::thread;

fn read_lines(file_path: &str) -> (Vec<String>, Vec<String>, Vec<String>) {
    let mut nodes: Vec<String> = vec![];
    let mut solids: Vec<String> = vec![];
    let mut shells: Vec<String> = vec![];
    let f = File::open(file_path);
    if let Err(e) = f {
        println!("{}", e);
        return (nodes, solids, shells);
    }
    let mut reader = BufReader::new(f.unwrap());
    let mut kwmode = KwBlock::Other;
    let mut line_buf = String::with_capacity(256);
    loop {
        line_buf.clear();
        let bytes_read = reader.read_line(&mut line_buf).unwrap();
        if bytes_read == 0 {
            break;
        }
        let line = line_buf.trim_end();
        if line.starts_with("*") {
            kwmode = if line.starts_with("*NODE") {
                KwBlock::Node
            } else if line.starts_with("*ELEMENT_SOLID") {
                KwBlock::ElementSolid
            } else if line.starts_with("*ELEMENT_SHELL") {
                KwBlock::ElementShell
            } else {
                KwBlock::Other
            };
            continue;
        }
        if line.starts_with("$") {
            continue;
        }
        match kwmode {
            KwBlock::Node => nodes.push(line.to_string()),
            KwBlock::ElementSolid => solids.push(line.to_string()),
            KwBlock::ElementShell => shells.push(line.to_string()),
            _ => (),
        }
    }
    (nodes, solids, shells)
}

pub fn read_k(
    file_path: &str,
) -> Result<
    (
        BTreeMap<usize, FENode>,
        BTreeMap<usize, FEElement>,
        BTreeMap<usize, FEElement>,
    ),
    Box<dyn std::error::Error>,
> {
    let (node_lines, solid_lines, shell_lines) = read_lines(file_path);
    // let nodes: BTreeMap<usize, LsNode> = BTreeMap::new();
    let nodes_handler = thread::spawn(move || {
        node_lines
            .par_iter()
            .map(|item| {
                let n: FENode = item.as_str().into();
                (n.n, n)
            })
            .collect()
    });
    let solids_handler = thread::spawn(move || {
        solid_lines
            .par_iter()
            .map(|item| {
                let mut e: FEElement = item.as_str().into();
                e.normalize_nodes();
                (e.n, e)
            })
            .collect()
    });
    let shells_handler = thread::spawn(move || {
        shell_lines
            .par_iter()
            .map(|item| {
                let mut e: FEElement = item.as_str().into();
                e.normalize_nodes();
                e.etype = ElementType::SHELL;
                (e.n, e)
            })
            .collect()
    });
    let nodes = nodes_handler.join().unwrap();
    let solids = solids_handler.join().unwrap();
    let shells = shells_handler.join().unwrap();
    Ok((nodes, solids, shells))
}
