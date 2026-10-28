use crate::utils::{Matrix3x3, determinant, mean_and_std, renumber_map};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::mem;
use std::path::PathBuf;
use std::str::FromStr;
use std::thread;

use crate::objects::feelement::ElementType;
use crate::objects::feelement::FEElement;
use crate::objects::fenode::FENode;

/// Скалярное поле по номеру узла
pub type NodeScalarField<'a> = (&'a str, &'a HashMap<usize, f32>);
/// Скалярное поле по номеру элемента (ключам `solids`/`shells`)
pub type CellScalarField<'a> = (&'a str, &'a HashMap<usize, f32>);
pub type SpacialHash = HashMap<(isize, isize, isize), Vec<usize>>;

/// Приводит имя поля к виду, который понимает legacy VTK:
/// пробелы заменяются на `_`, пустое имя становится `field_<номер>`.
fn vtk_field_name(name: &str) -> String {
    let name = name.trim().replace([' ', '\t'], "_");
    if name.is_empty() {
        "field".to_string()
    } else {
        name
    }
}

enum KwBlock {
    Node,
    ElementSolid,
    ElementShell,
    Other,
}

// /// Тип узла сетки
// #[derive(Debug, Clone, Copy, PartialEq)]
// pub enum BoundaryKind {
//     Interior,
//     Face,   // соседи лежат в полупространстве (грань)
//     Edge,   // соседи лежат в двугранном угле (ребро)
//     Corner, // соседи лежат в трёхгранном угле (угол)
// }

fn is_contains(st: &HashSet<usize>, val: &usize) -> bool {
    if st.is_empty() {
        true
    } else {
        st.contains(&val)
    }
}

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

pub struct FEModel {
    pub nodes: BTreeMap<usize, FENode>,
    pub solids: BTreeMap<usize, FEElement>,
    pub shells: BTreeMap<usize, FEElement>,
    pub tolerance: f64,
}

impl Default for FEModel {
    fn default() -> Self {
        Self {
            nodes: BTreeMap::new(),
            solids: BTreeMap::new(),
            shells: BTreeMap::new(),
            tolerance: 1e-6,
        }
    }
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
                        "c3d8" | "c3d8r" => KwBlock::ElementSolid,
                        "ax4" => KwBlock::ElementShell,
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

impl FEModel {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read(&mut self, path: &str) -> Result<(), String> {
        let path = PathBuf::from_str(path).map_err(|e| e.to_string())?;
        if !path.exists() {
            return Err("File not found".to_string());
        }
        let extension = path
            .extension()
            .map(|ext| ext.to_string_lossy().to_lowercase())
            .unwrap_or_default();
        let file_path = path.to_string_lossy().to_string();
        match extension.as_str() {
            "k" => {
                let (nodes, solids, shells) = read_k(&file_path).map_err(|e| e.to_string())?;
                self.nodes = nodes;
                self.solids = solids;
                self.shells = shells;
            }
            "inp" => {
                let (nodes, solids, shells) =
                    read_flat_abaqus_input(&file_path).map_err(|e| e.to_string())?;
                self.nodes = nodes;
                self.solids = solids;
                self.shells = shells;
            }
            other => {
                return Err(format!(
                    "Unsupported file extension \"{other}\", expected \"k\" or \"inp\""
                ));
            }
        }
        Ok(())
    }

    /// Координаты центра элемента с номером n
    pub fn element_center(&self, en: &usize) -> Option<[f64; 3]> {
        let el = self.solids.get(en)?;
        let n = el.nodes_count();
        let mut crds = (0..n).into_iter().fold([0.0, 0.0, 0.0], |acc, n| {
            let x = acc[0] + self.nodes[&el.nodes[n]].x;
            let y = acc[1] + self.nodes[&el.nodes[n]].y;
            let z = acc[2] + self.nodes[&el.nodes[n]].z;
            [x, y, z]
        });
        crds[0] /= n as f64;
        crds[1] /= n as f64;
        crds[2] /= n as f64;
        Some(crds)
    }

    pub fn element_volume(&self, el_number: usize) -> Option<f64> {
        if let Some(el) = self.solids.get(&el_number) {
            if el.etype != ElementType::SOLID {
                return Some(0f64);
            }
            match el.nodes_count() {
                4 => {
                    let mut idxs = [0usize; 4];
                    for i in 0..4 {
                        idxs[i] = el.nodes[FOUR_NODES_INDEXES[0][i] - 1];
                    }
                    let v = self.vol4nodes(&idxs)?;
                    Some(v)
                }
                6 => {
                    let mut v = 0.0;
                    let mut idxs = [0usize; 4];
                    for idx_row in SIX_NODES_INDEXES {
                        for j in 0..4 {
                            idxs[j] = el.nodes[idx_row[j] - 1];
                        }
                        let vv = self.vol4nodes(&idxs)?;
                        v += vv;
                    }
                    Some(v)
                }
                8 => {
                    let mut v = 0.0;
                    let mut idxs = [0usize; 4];
                    for idx_row in EIGHT_NODES_INDEXES {
                        for j in 0..4 {
                            idxs[j] = el.nodes[idx_row[j] - 1];
                        }
                        let vv = self.vol4nodes(&idxs)?;
                        v += vv;
                    }
                    Some(v)
                }
                _ => None,
            }
        } else {
            None
        }
    }

    fn vol4nodes(&self, nodes_numbers: &[usize; 4]) -> Option<f64> {
        let n1 = self.nodes.get(&nodes_numbers[0])?.coordinates();
        let n2 = self.nodes.get(&nodes_numbers[1])?.coordinates();
        let n3 = self.nodes.get(&nodes_numbers[2])?.coordinates();
        let n4 = self.nodes.get(&nodes_numbers[3])?.coordinates();
        let mut matr: Matrix3x3 = Matrix3x3::default();
        for i in 0..3 {
            matr[0][i] = n2[i] - n1[i];
        }
        for i in 0..3 {
            matr[1][i] = n3[i] - n1[i];
        }
        for i in 0..3 {
            matr[2][i] = n4[i] - n1[i];
        }
        Some(determinant(&matr).abs() / 6.0)
    }
    pub fn calculate_element_volumes(&mut self) {
        let volumes = self
            .solids
            .par_iter()
            .map(|(k, _)| (*k, self.element_volume(*k).unwrap_or_default()))
            .collect::<BTreeMap<usize, f64>>();
        for (k, v) in volumes {
            self.solids.get_mut(&k).unwrap().volume = v;
        }
    }
    pub fn elements_by_part(&self, part: &usize) -> Vec<usize> {
        self.solids
            .par_iter()
            .filter_map(|(k, e)| if e.part == *part { Some(*k) } else { None })
            .collect()
    }

    /// Method for save model to k-file
    ///
    /// *file_name* is the path to k-file
    pub fn save(&mut self, file_name: &str) -> Result<(), std::io::Error> {
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(f);
        writer.write_all(b"*KEYWORD\n")?;
        if !self.nodes.is_empty() {
            writer.write_all(b"*NODE\n")?;
            let nodes = self
                .nodes
                .par_iter()
                .map(|(_, n)| n.into())
                .collect::<Vec<String>>()
                .join("\n");
            writeln!(writer, "{}", nodes)?;
        }
        if !self.solids.is_empty() {
            writer.write_all(b"*ELEMENT_SOLID\n")?;
            let solids = self
                .solids
                .par_iter()
                .map(|(_, e)| e.try_into().unwrap())
                .collect::<Vec<String>>()
                .join("\n");
            writeln!(writer, "{}", solids)?;
        }
        if !self.shells.is_empty() {
            writer.write_all(b"*ELEMENT_SHELL\n")?;
            let shells = self
                .shells
                .par_iter()
                .map(|(_, e)| e.try_into().unwrap())
                .collect::<Vec<String>>()
                .join("\n");
            writeln!(writer, "{}", shells)?;
        }
        writer.write_all(b"*END")?;
        Ok(())
    }

    /// Сохраняет модель в legacy VTK (ASCII) со скалярными полями в узлах и элементах.
    ///
    /// *point_fields* — узловые поля: (имя, &значения по номеру узла).
    /// *cell_fields*  — поля по ячейкам: (имя, &значения по номеру элемента, т.е. по ключу `solids`/`shells`).
    /// Для отсутствующих узлов/элементов пишется `NaN`. Порядок значений в файле — порядок
    /// точек (`nodes`) и ячеек (`solids`, затем `shells`), как в блоках `POINTS`/`CELLS`.
    ///
    /// mesh.save_vtk("mesh.vtk", None, None); запись без полей
    /// mesh.save_vtk("mesh.vtk", Some(&[("x", &x), ("y", &y)]), None); только узловые поля
    /// mesh.save_vtk("mesh.vtk", None, Some(&[("volume", &volumes)])); только по ячейкам
    /// mesh.save_vtk("mesh.vtk", Some(&[("x", &x)]), Some(&[("volume", &volumes)])); и то, и другое
    ///
    /// При заданных полях дополнительно пишутся служебные ячейковые поля
    /// `element_index` и `part` (сопоставление ячейки VTK с элементом модели).
    pub fn save_vtk(
        &self,
        file_name: &str,
        point_fields: Option<&[NodeScalarField]>,
        cell_fields: Option<&[CellScalarField]>,
    ) -> Result<(), std::io::Error> {
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(f);

        // 1. Сплошная нумерация узлов 0..n, номер узла -> индекс точки
        let index: HashMap<usize, usize> = self
            .nodes
            .keys()
            .enumerate()
            .map(|(i, n)| (*n, i))
            .collect();

        // 2. Элементы -> (vtk_type, Vec<индексы точек>, номер элемента, part)
        let mut cells: Vec<(u8, Vec<usize>, usize, usize)> = Vec::new();
        for e in self.solids.values().chain(self.shells.values()) {
            let ids: Vec<usize> = e.nodes[..e.nodes_count()]
                .iter()
                .filter_map(|n| index.get(n).copied())
                .collect();
            let vtk_type = match (e.etype.clone(), ids.len()) {
                (ElementType::SOLID, 4) => 10, // VTK_TETRA
                (ElementType::SOLID, 6) => 13, // VTK_WEDGE
                (ElementType::SOLID, 8) => 12, // VTK_HEXAHEDRON
                (ElementType::SHELL, 3) => 5,  // VTK_TRIANGLE
                (ElementType::SHELL, 4) => 9,  // VTK_QUAD
                _ => continue,
            };
            cells.push((vtk_type, ids, e.n, e.part));
        }
        let total_conn: usize = cells.iter().map(|(_, ids, _, _)| ids.len()).sum();

        // 3. Заголовок
        writeln!(writer, "# vtk DataFile Version 3.0")?;
        writeln!(writer, "lsmesh export")?;
        writeln!(writer, "ASCII")?;
        writeln!(writer, "DATASET UNSTRUCTURED_GRID")?;

        // 4. Точки
        writeln!(writer, "POINTS {} double", self.nodes.len())?;
        for n in self.nodes.values() {
            writeln!(writer, "{} {} {}", n.x, n.y, n.z)?;
        }

        // 5. Ячейки: каждая строка "k id0 id1 ... id(k-1)"
        writeln!(writer, "CELLS {} {}", cells.len(), cells.len() + total_conn)?;
        for (_, ids, _, _) in &cells {
            write!(writer, "{}", ids.len())?;
            for id in ids {
                write!(writer, " {}", id)?;
            }
            writeln!(writer)?;
        }
        writeln!(writer, "CELL_TYPES {}", cells.len())?;
        for (t, _, _, _) in &cells {
            writeln!(writer, "{}", t)?;
        }

        // 6. Узловые скалярные поля (по порядку вывода точек)
        if let Some(fields) = point_fields {
            if !fields.is_empty() && !self.nodes.is_empty() {
                writeln!(writer, "POINT_DATA {}", self.nodes.len())?;
                for (name, values) in fields {
                    writeln!(writer, "SCALARS {} double 1", vtk_field_name(name))?;
                    writeln!(writer, "LOOKUP_TABLE default")?;
                    for n in self.nodes.keys() {
                        writeln!(writer, "{}", values.get(n).copied().unwrap_or(f32::NAN))?;
                    }
                }
            }
        }

        // 7. Служебные поля - номера ячеек и номера частей
        if !cells.is_empty() {
            writeln!(writer, "CELL_DATA {}", cells.len())?;
            writeln!(writer, "SCALARS element_index int 1")?;
            writeln!(writer, "LOOKUP_TABLE default")?;
            for (_, _, n, _) in &cells {
                writeln!(writer, "{}", n)?;
            }
            writeln!(writer, "SCALARS part int 1")?;
            writeln!(writer, "LOOKUP_TABLE default")?;
            for (_, _, _, p) in &cells {
                writeln!(writer, "{}", p)?;
            }
        }

        // 8. Скалярные поля по ячейкам (по порядку вывода ячеек)
        if let Some(fields) = cell_fields {
            if !fields.is_empty() && !cells.is_empty() {
                for (name, values) in fields {
                    writeln!(writer, "SCALARS {} double 1", vtk_field_name(name))?;
                    writeln!(writer, "LOOKUP_TABLE default")?;
                    for (_, _, n, _) in &cells {
                        writeln!(writer, "{}", values.get(n).copied().unwrap_or(f32::NAN))?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Translates model (all nodes) by vector v
    pub fn translate(&mut self, v: [f64; 3]) {
        self.nodes.iter_mut().for_each(|(_, n)| n.translate(&v));
    }

    /// Determinate BBox of model
    pub fn bbox(&self) -> BBox {
        if self.nodes.is_empty() {
            return BBox::default();
        }
        self.nodes
            .values()
            .map(|n| BBox {
                min_x: n.x,
                max_x: n.x,
                min_y: n.y,
                max_y: n.y,
                min_z: n.z,
                max_z: n.z,
            })
            .reduce(|a, b| BBox {
                min_x: a.min_x.min(b.min_x),
                max_x: a.max_x.max(b.max_x),
                min_y: a.min_y.min(b.min_y),
                max_y: a.max_y.max(b.max_y),
                min_z: a.min_z.min(b.min_z),
                max_z: a.max_z.max(b.max_z),
            })
            .unwrap_or_default()
    }

    /// Volume of model bbox
    pub fn bbox_volume(&self) -> f64 {
        let bb = self.bbox();
        (bb.max_x - bb.min_x) * (bb.max_y - bb.min_y) * (bb.max_z - bb.min_z)
    }

    pub fn set_tolerance(&mut self, tolerance: f64) {
        self.tolerance = tolerance;
    }

    pub fn renumber_solids(&mut self, start_number: usize) {
        if !self.solids.is_empty() {
            let solids = mem::take(&mut self.solids);
            (self.solids, _) = renumber_map(solids, start_number, false);
        }
    }

    pub fn renumber_shells(&mut self, start_number: usize) {
        if !self.shells.is_empty() {
            let shells = mem::take(&mut self.shells);
            (self.shells, _) = renumber_map(shells, start_number, false);
        }
    }

    pub fn renumber_elements(&mut self, start_number: usize) {
        if !self.solids.is_empty() {
            self.renumber_solids(start_number);
        }
        if !self.shells.is_empty() {
            self.renumber_shells(start_number + self.solids.len());
        }
    }

    pub fn renumber_nodes(&mut self, start_number: usize) {
        let nodes = mem::take(&mut self.nodes);
        let mut new_numbers: HashMap<usize, usize>;
        (self.nodes, new_numbers) = renumber_map(nodes, start_number, true);
        new_numbers.insert(0, 0);
        if !self.solids.is_empty() {
            self.solids.iter_mut().for_each(|(_k, s)| {
                for n in s.nodes.iter_mut() {
                    *n = new_numbers[n];
                }
            })
        }
        if !self.shells.is_empty() {
            self.shells.iter_mut().for_each(|(_k, s)| {
                for n in s.nodes.iter_mut() {
                    *n = new_numbers[n];
                }
            })
        }
    }
    pub fn parts(&self) -> Vec<usize> {
        let mut parts: HashSet<usize> = HashSet::new();
        self.solids.iter().for_each(|(_k, s)| {
            parts.insert(s.part);
        });
        self.shells.iter().for_each(|(_k, s)| {
            parts.insert(s.part);
        });
        let mut parts: Vec<usize> = parts.into_iter().collect();
        parts.sort();
        parts
    }

    pub fn renumber_parts(&mut self, start_number: usize) {
        let new_numbers: HashMap<usize, usize> = self
            .parts()
            .iter()
            .enumerate()
            .map(|(i, p)| (*p, i + start_number))
            .collect();
        self.solids.iter_mut().for_each(|(_k, v)| {
            v.part = new_numbers[&v.part];
        });
        self.shells.iter_mut().for_each(|(_k, v)| {
            v.part = new_numbers[&v.part];
        });
    }

    /// Searchs for unconnected parts of FE model (groups of elements whith no common nodes to elements in other group)
    /// if parts = &[] searchs in all parts of the model
    ///    parts = &[1] searchs only in part 1
    /// returns Vec<Vec<usize>> - vector of vector of element indexes
    /// [[1,2], [3,4]] - elements 1 and 2 in first region, 3 and 4 in second
    pub fn find_unconnected_regions(&self, parts: &[usize]) -> Vec<Vec<usize>> {
        let parts: HashSet<usize> = parts.iter().copied().collect();
        let mut nodes_to_elements = HashMap::<usize, Vec<usize>>::new();
        for (eid, e) in self
            .solids
            .iter()
            .filter(|m| is_contains(&parts, &m.1.part))
        {
            for nid in e.nodes {
                if nid == 0 {
                    break;
                }
                nodes_to_elements
                    .entry(nid)
                    .or_insert_with(Vec::new)
                    .push(*eid);
            }
        }
        let mut visited = HashSet::<usize>::new();
        let mut components = Vec::<Vec<usize>>::new();
        for (start, _) in self
            .solids
            .iter()
            .filter(|m| is_contains(&parts, &m.1.part))
        {
            if visited.contains(start) {
                continue;
            }
            let mut stack = vec![start];
            visited.insert(*start);
            let mut comp = Vec::<usize>::new();
            while stack.len() > 0 {
                let e = stack.pop().unwrap();
                comp.push(*e);
                for nid in self.solids[e].nodes {
                    if nid == 0 {
                        break;
                    }
                    for ee in nodes_to_elements[&nid].iter() {
                        if !visited.contains(&ee) {
                            visited.insert(*ee);
                            stack.push(&ee);
                        }
                    }
                }
            }
            components.push(comp);
        }
        components
    }

    /// возвращает пространственный hash модели в виде словаря
    /// ключи - индексы клеток размером cell_size
    /// значение - вектор индексов узлов, попавших в клетку
    pub fn nodes_spacial_hash(&self, cell_size: f32) -> SpacialHash {
        let mut hash: SpacialHash = HashMap::new();
        for n in self.nodes.values() {
            let node_hash = n.spacial_hash(cell_size);
            hash.entry(node_hash).or_insert_with(Vec::new).push(n.n);
        }
        hash
    }

    /// Ищет соседние узлы внутри сферы радиуса radius
    /// Возвращет индексы узлов
    /// Пространственный хэш *spacial_hash* должен быть построен для того же размера ячейки *cell_size*,
    /// что передается в функцию
    pub fn node_neibors(
        &self,
        n0: &FENode,
        spacial_hash: &SpacialHash,
        cell_size: f32,
        radius: f32,
    ) -> Vec<usize> {
        let mut result: Vec<usize> = Vec::with_capacity(64);
        let n_cells = (radius / cell_size).ceil() as isize;
        let n0_sh = n0.spacial_hash(cell_size);
        let r2 = radius * radius;
        for i in (n0_sh.0 - n_cells)..=(n0_sh.0 + n_cells) {
            for j in (n0_sh.1 - n_cells)..=(n0_sh.1 + n_cells) {
                for k in (n0_sh.2 - n_cells)..=(n0_sh.2 + n_cells) {
                    if let Some(cell_nodes) = spacial_hash.get(&(i, j, k)) {
                        for &n in cell_nodes {
                            let n1 = &self.nodes[&n];
                            let dx = n0.x - n1.x;
                            let dy = n0.y - n1.y;
                            let dz = n0.z - n1.z;
                            if (dx * dx + dy * dy + dz * dz) as f32 <= r2 {
                                result.push(n);
                            }
                        }
                    }
                }
            }
        }
        result
    }

    pub fn node_neibors_ellipse(
        &self,
        n0: &FENode,
        spacial_hash: &SpacialHash,
        cell_size: f32,
        rx: f32,
        ry: f32,
        rz: f32,
    ) -> Vec<usize> {
        let mut result: Vec<usize> = Vec::with_capacity(64);
        let n_cells_x = (rx / cell_size).ceil() as isize;
        let n_cells_y = (ry / cell_size).ceil() as isize;
        let n_cells_z = (rz / cell_size).ceil() as isize;
        let n0_sh = n0.spacial_hash(cell_size);
        let rx2 = (rx * rx) as f64;
        let ry2 = (ry * ry) as f64;
        let rz2 = (rz * rz) as f64;
        for i in (n0_sh.0 - n_cells_x)..=(n0_sh.0 + n_cells_x) {
            for j in (n0_sh.1 - n_cells_y)..=(n0_sh.1 + n_cells_y) {
                for k in (n0_sh.2 - n_cells_z)..=(n0_sh.2 + n_cells_z) {
                    if let Some(cell_nodes) = spacial_hash.get(&(i, j, k)) {
                        for &n in cell_nodes {
                            let n1 = &self.nodes[&n];
                            let dx = n0.x - n1.x;
                            let dy = n0.y - n1.y;
                            let dz = n0.z - n1.z;
                            if (dx * dx / rx2 + dy * dy / ry2 + dz * dz / rz2) as f32 <= 1.0 {
                                result.push(n);
                            }
                        }
                    }
                }
            }
        }
        result
    }

    pub fn node_neibors_with_distances2(
        &self,
        n0: &FENode,
        spacial_hash: &SpacialHash,
        cell_size: f32,
        radius: f32,
    ) -> HashMap<usize, f32> {
        let mut result: HashMap<usize, f32> = HashMap::with_capacity(64);
        let n_cells = (radius / cell_size).ceil() as isize;
        let n0_sh = n0.spacial_hash(cell_size);
        let r2 = radius * radius;
        for i in (n0_sh.0 - n_cells)..=(n0_sh.0 + n_cells) {
            for j in (n0_sh.1 - n_cells)..=(n0_sh.1 + n_cells) {
                for k in (n0_sh.2 - n_cells)..=(n0_sh.2 + n_cells) {
                    if let Some(cell_nodes) = spacial_hash.get(&(i, j, k)) {
                        for &n in cell_nodes {
                            let n1 = &self.nodes[&n];
                            let dx = n0.x - n1.x;
                            let dy = n0.y - n1.y;
                            let dz = n0.z - n1.z;
                            let d2 = (dx * dx + dy * dy + dz * dz) as f32;
                            if d2 <= r2 {
                                result.insert(n, d2);
                            }
                        }
                    }
                }
            }
        }
        result
    }
}

impl fmt::Display for FEModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Number of nodes: {}, number of solids: {}, number of shells: {}",
            self.nodes.len(),
            self.solids.len(),
            self.shells.len()
        )
    }
}

#[derive(Debug, Default)]
pub struct BBox {
    pub min_x: f64,
    pub max_x: f64,
    pub min_y: f64,
    pub max_y: f64,
    pub min_z: f64,
    pub max_z: f64,
}

const FOUR_NODES_INDEXES: [[usize; 4]; 1] = [[1, 2, 3, 4]];
const SIX_NODES_INDEXES: [[usize; 4]; 3] = [[1, 2, 4, 5], [5, 2, 4, 3], [5, 4, 3, 6]];
const EIGHT_NODES_INDEXES: [[usize; 4]; 5] = [
    [1, 5, 6, 8],
    [3, 7, 6, 8],
    [1, 3, 4, 8],
    [1, 3, 6, 8],
    [1, 2, 3, 6],
];
