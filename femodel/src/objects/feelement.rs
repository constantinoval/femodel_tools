use crate::utils::{HasNum, parse_sub_str};
use std::convert::TryInto;
use std::fmt::Write;

#[derive(Debug, Default, PartialEq, Clone)]
pub enum ElementType {
    #[default]
    SOLID,
    SHELL,
    BEAM,
    SPH,
}

#[derive(Debug, Default, Clone)]
pub struct FEElement {
    pub n: usize,
    pub part: usize,
    pub nodes: [usize; 8],
    pub volume: f64,
    pub etype: ElementType,
}

impl From<&str> for FEElement {
    fn from(value: &str) -> Self {
        if value.as_bytes().contains(&b',') {
            let v: Vec<_> = value.split(',').map(|substr| substr.trim()).collect();
            let n = v
                .first()
                .unwrap_or(&"")
                .parse::<usize>()
                .unwrap_or_default();
            let part = v.get(1).unwrap_or(&"").parse::<usize>().unwrap_or_default();
            let nodes: [usize; 8] = (2..10)
                .map(|i| v.get(i).unwrap_or(&"").parse::<usize>().unwrap_or_default())
                .collect::<Vec<usize>>()
                .try_into()
                .unwrap();
            FEElement {
                n,
                part,
                nodes,
                volume: 0.0f64,
                etype: ElementType::SOLID,
            }
        } else {
            let n: usize = parse_sub_str(value, 0, 8);
            let part: usize = parse_sub_str(value, 8, 16);
            let nodes: [usize; 8] = (0..8)
                .map(|i| parse_sub_str::<usize>(value, 16 + i * 8, 16 + (i + 1) * 8))
                .collect::<Vec<usize>>()
                .try_into()
                .unwrap();
            FEElement {
                n,
                part,
                nodes,
                volume: 0.0f64,
                etype: ElementType::SOLID,
            }
        }
    }
}

impl FEElement {
    pub fn normalize_nodes(&mut self) {
        let mut normalized = [0usize; 8];
        let mut position = 0;
        'l: for i in 0..self.nodes.len() {
            for j in 0..i {
                if self.nodes[i] == self.nodes[j] {
                    continue 'l;
                }
            }
            normalized[position] = self.nodes[i];
            position += 1;
        }
        self.nodes = normalized;
    }

    pub fn nodes_count(&self) -> usize {
        for i in 0..self.nodes.len() {
            if self.nodes[i] == 0 {
                return i;
            }
        }
        self.nodes.len()
    }
}

impl TryFrom<&FEElement> for String {
    type Error = String;
    fn try_from(el: &FEElement) -> Result<String, String> {
        let mut s: String = String::with_capacity(80);
        write!(&mut s, " {:7} {:7}", el.n, el.part).unwrap();
        match el.etype {
            ElementType::SOLID => match el.nodes_count() {
                8 => {
                    write!(
                        &mut s,
                        "{}",
                        &el.nodes
                            .iter()
                            .map(|n| format!(" {:7}", n))
                            .collect::<Vec<String>>()
                            .join(""),
                    )
                    .unwrap();
                    Ok(s)
                }
                4 => {
                    write!(
                        &mut s,
                        "{}",
                        &el.nodes[..4]
                            .iter()
                            .map(|n| format!(" {:7}", n))
                            .collect::<Vec<String>>()
                            .join(""),
                    )
                    .unwrap();
                    for _ in 0..4 {
                        write!(&mut s, " {:7}", el.nodes[3]).unwrap();
                    }
                    Ok(s)
                }
                6 => {
                    write!(
                        &mut s,
                        "{}",
                        &el.nodes[..4]
                            .iter()
                            .map(|n| format!(" {:7}", n))
                            .collect::<Vec<String>>()
                            .join(""),
                    )
                    .unwrap();
                    write!(&mut s, " {0:7} {0:7} {1:7} {1:7}", el.nodes[4], el.nodes[5]).unwrap();
                    Ok(s)
                }
                _ => Err("Wrong number of nodes!".to_string()),
            },
            ElementType::SHELL => match el.nodes_count() {
                4 => {
                    write!(
                        &mut s,
                        "{}",
                        &el.nodes[..4]
                            .iter()
                            .map(|n| format!(" {:7}", n))
                            .collect::<Vec<String>>()
                            .join(""),
                    )
                    .unwrap();
                    Ok(s)
                }
                3 => {
                    write!(
                        &mut s,
                        "{}",
                        &el.nodes[..3]
                            .iter()
                            .map(|n| format!(" {:7}", n))
                            .collect::<Vec<String>>()
                            .join(""),
                    )
                    .unwrap();
                    write!(&mut s, " {:7}", el.nodes[2]).unwrap();
                    Ok(s)
                }
                _ => Err("Wrong number of nodes!".to_string()),
            },
            _ => Err("Not imlemented element type".to_string()),
        }
    }
}

impl HasNum for FEElement {
    fn mut_num(&mut self) -> &mut usize {
        &mut self.n
    }
}
