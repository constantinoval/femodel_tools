use crate::objects::feelement::ElementType;
use crate::prelude::*;
use crate::readers::common_types::*;
use anyhow::Context;
use rayon::prelude::*;
use regex::Regex;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug)]
enum ElsetType {
    List(Vec<usize>),
    Generate((usize, usize, usize)),
}

#[derive(Debug)]
struct Elset {
    name: String,
    set_type: ElsetType,
}

impl Elset {
    fn from_strings(lines: Vec<String>) -> Option<Self> {
        let elset_r =
            Regex::new(r"(?i:^*elset, (elset=(?P<name>[^\s,]*))?[,\s]*(?P<generate>generate)?).*")
                .unwrap();
        let c = elset_r.captures(lines.first()?)?;
        let name = c.name("name")?.as_str();
        if let Some(_generate) = c.name("generate") {
            let params: Vec<usize> = lines
                .get(1)?
                .split(',')
                .filter_map(|s| s.trim().parse::<usize>().ok())
                .collect();
            if params.len() == 3 {
                Some(Elset {
                    name: name.to_string(),
                    set_type: ElsetType::Generate((params[0], params[1], params[2])),
                })
            } else {
                None
            }
        } else {
            let list = lines[1..]
                .iter()
                .flat_map(|line| line.split(','))
                .filter_map(|s| s.trim().parse::<usize>().ok())
                .collect();
            Some(Elset {
                name: name.to_string(),
                set_type: ElsetType::List(list),
            })
        }
    }

    fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        match &self.set_type {
            ElsetType::List(v) => {
                Box::new(v.iter().copied()) as Box<dyn Iterator<Item = usize> + '_>
            }
            ElsetType::Generate((start, stop, step)) => {
                Box::new((*start..=*stop).step_by(*step.max(&1)))
            }
        }
    }
}

fn get_section_elset_name(line: &str) -> Option<String> {
    let section_r = Regex::new(r"(?i:^*solid section, (elset=(?P<elset>[^\s,]*))?.*)").unwrap();
    Some(
        section_r
            .captures(line)?
            .name("elset")?
            .as_str()
            .to_string(),
    )
}

/// Чтение сетки из файла Abaqus (формат без parts и assembly)
pub fn read_flat_abaqus_input(
    file_path: &str,
) -> anyhow::Result<(
    BTreeMap<usize, FENode>,
    BTreeMap<usize, FEElement>,
    BTreeMap<usize, FEElement>,
)> {
    fn split_line(line: &str, sep: char) -> Vec<&str> {
        line.split(sep)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect()
    }

    fn parse_node(line: &str) -> Option<FENode> {
        let parts = split_line(line, ',');
        let n = parts.first()?.parse::<usize>().ok()?;
        if n == 0 {
            return None;
        };
        let x: f64 = parts
            .get(1)
            .unwrap_or(&"")
            .parse::<f64>()
            .unwrap_or_default();
        let y: f64 = parts
            .get(2)
            .unwrap_or(&"")
            .parse::<f64>()
            .unwrap_or_default();
        let z: f64 = parts
            .get(3)
            .unwrap_or(&"")
            .parse::<f64>()
            .unwrap_or_default();
        Some(FENode { n, x, y, z })
    }

    fn parse_element(line: &str, etype: ElementType) -> Option<FEElement> {
        let parts = split_line(line, ',');
        let n = parts.first()?.parse::<usize>().ok()?;
        if n == 0 {
            return None;
        }
        let v: Vec<usize> = parts[1..]
            .iter()
            .map(|s| s.parse::<usize>().unwrap_or(0))
            .collect();
        let mut nodes = [0usize; 8];
        let nc = v.len().min(8);
        nodes[..nc].copy_from_slice(&v[..nc]);
        Some(FEElement {
            n,
            part: 1,
            etype,
            nodes,
            volume: 0.0,
        })
    }

    let mut node_lines: Vec<String> = Vec::with_capacity(100000);
    let mut solid_lines: Vec<String> = Vec::with_capacity(100000);
    let mut shell_lines: Vec<String> = Vec::with_capacity(100000);

    let mut bt = KwBlock::Other;
    let f = File::open(file_path)?;
    let reader = BufReader::new(&f);
    for l in reader.lines() {
        let line = l.unwrap_or_default();
        if line.starts_with("**") {
            continue;
        }
        if line.starts_with("*") {
            let kw = line[1..]
                .split(',')
                .map(|s| s.trim().to_lowercase())
                .collect::<Vec<String>>();
            bt = match kw[0].as_str() {
                "node" => KwBlock::Node,
                "element" => {
                    let etype = kw[1].split("=").last().unwrap_or("").trim();
                    match etype {
                        "c3d8" | "c3d8r" | "c3d4" => KwBlock::ElementSolid,
                        "cax4r" => KwBlock::ElementShell,
                        _ => KwBlock::Other,
                    }
                }
                _ => KwBlock::Other,
            };
            continue;
        }
        match bt {
            KwBlock::Node => node_lines.push(line),
            KwBlock::ElementSolid => solid_lines.push(line),
            KwBlock::ElementShell => shell_lines.push(line),
            KwBlock::Other => {}
        }
    }

    let nodes: BTreeMap<usize, FENode> = node_lines
        .par_iter()
        .filter_map(|l| {
            let nd = parse_node(l)?;
            Some((nd.n, nd))
        })
        .collect();

    let solids: BTreeMap<usize, FEElement> = solid_lines
        .par_iter()
        .filter_map(|l| {
            let e = parse_element(l, ElementType::SOLID)?;
            Some((e.n, e))
        })
        .collect();

    let shells: BTreeMap<usize, FEElement> = shell_lines
        .par_iter()
        .filter_map(|l| {
            let e = parse_element(l, ElementType::SHELL)?;
            Some((e.n, e))
        })
        .collect();
    Ok((nodes, solids, shells))
}

/// Читает из inp файла скалярные поля записанные ключевым словом *Initial conditions, type=FIELD
pub fn read_nodal_fields_from_inp(
    file_path: &str,
) -> anyhow::Result<HashMap<usize, HashMap<usize, f32>>> {
    fn process_data_lines(data_lines: &Vec<String>) -> HashMap<usize, f32> {
        data_lines
            .par_iter()
            .filter_map(|line| {
                let mut it = line.split(',');
                let n = it.next()?.trim().parse::<usize>().ok()?;
                let f = it.next()?.trim().parse::<f32>().ok()?;
                Some((n, f))
            })
            .collect()
    }

    let mut result: HashMap<usize, HashMap<usize, f32>> = HashMap::new();
    let f = File::open(file_path)?;
    let reader = BufReader::new(f);
    let mut inside_block = false;
    let mut data_lines: Vec<String> = Vec::new();
    let mut field_id: usize = 0;
    let ic_regex =
        Regex::new(r"^*Initial conditions,\s*type=FIELD,\s*variable\s*=\s*(?P<id>\d+).*").unwrap();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("**") {
            continue;
        }
        if line.starts_with("*") {
            if !data_lines.is_empty() {
                result.insert(field_id, process_data_lines(&data_lines));
                data_lines.clear();
            }
            if let Some(id) = ic_regex
                .captures(&line)
                .context("Can't capture")?
                .name("id")
            {
                field_id = id.as_str().parse::<usize>()?;
                data_lines.clear();
                inside_block = true;
            } else {
                inside_block = false;
            }
        }
        if inside_block {
            data_lines.push(line);
        }
    }
    if !data_lines.is_empty() {
        result.insert(field_id, process_data_lines(&data_lines));
        data_lines.clear();
    }
    Ok(result)
}
