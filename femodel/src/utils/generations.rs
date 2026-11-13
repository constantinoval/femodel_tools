use crate::objects::femodel::FEModel;
use crate::utils::mean_and_std;
use anyhow::Context;
use fast_math;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rand_distr::{Distribution, StandardNormal};
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::hash_map::Keys;
use std::fmt;
use std::io::Write;

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug)]
pub enum FieldType {
    Isotropic(f32),
    Orthotropic(f32, f32, f32),
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug)]
pub struct FieldProperties {
    pub mean: f32,
    pub std: f32,
    pub id: u16,
    pub name: String,
    pub field_type: FieldType,
    pub inverse: bool,
}

// impl FieldProperties {
//     pub fn new() -> Self {
//         Self {
//             mean: 0.0,
//             std: 1.0,
//             id: 1,
//             name: "".to_string(),
//             field_type: FieldType::Isotropic(1.0),
//             inverse: false,
//         }
//     }
// }

// pub struct GeneratedField {
//     pub field_property: FieldProperties,
//     pub field: HashMap<usize, f32>,
// }

impl From<&FieldType> for String {
    fn from(value: &FieldType) -> Self {
        match value {
            FieldType::Isotropic(l) => format!("isotropic(Lc={})", l).to_string(),
            FieldType::Orthotropic(lx, ly, lz) => {
                format!("orthotropic(Lx={}, Ly={}, Lz={})", lx, ly, lz).to_string()
            }
        }
    }
}

// impl From<&FieldProperties> for String {
//     fn from(value: &FieldProperties) -> Self {
//         let mut result = "".to_string();
//         result.push_str(
//             &format!(
//                 "Field number {}, mean={}, std={}, correlation length: {}",
//                 &value.id,
//                 &value.mean,
//                 &value.std,
//                 String::from(&value.field_type)
//             )
//             .to_string(),
//         );
//         result
//     }
// }

impl FEModel {
    /// Генерирует изотропное коррелированное узловое случайное поле с длиной корреляции 2*l
    pub fn generate_nodal_correlated_random_field(&self, l: f32) -> HashMap<usize, f32> {
        let bar = ProgressBar::new(self.nodes.len() as u64).with_style(
            ProgressStyle::with_template(
                "{prefix} [{bar:40.cyan/blue}] {pos:>9}/{len:9} {percent:>3}% ETA {eta:>5}",
            )
            .unwrap()
            .progress_chars("=> "),
        );
        let r_cut = 3.0 * l;
        let cell_size = r_cut;
        let f: HashMap<usize, f32> = self
            .nodes
            .par_iter()
            .map(|(&n, _)| {
                let mut rng = rand::rng();
                (n, StandardNormal.sample(&mut rng))
            })
            .collect();
        let sh = self.nodes_spacial_hash(cell_size);
        let z: HashMap<usize, f32> = self
            .nodes
            .par_iter()
            .progress_with(bar.clone())
            .map(|(_, n0)| {
                let nds = self.node_neibors_with_distances2(n0, &sh, cell_size, r_cut);
                let (s1, s2, s2q) =
                    nds.iter()
                        .fold((0.0f32, 0.0f32, 0.0f32), |(a1, a2, a3), (&m, d)| {
                            let w = fast_math::exp(-d / 2.0 / l / l);
                            (a1 + w * f[&m], a2 + w, a3 + w * w)
                        });
                let mean = s1 / s2;
                let sigma = (s2q / (s2 * s2)).sqrt();
                (n0.n, mean / sigma)
            })
            .collect();
        let zs = z.values().copied().collect::<Vec<f32>>();
        let (mean, std) = mean_and_std(&zs);
        let z: HashMap<usize, f32> = z
            .par_iter()
            .map(|(n, zz)| {
                let zzz = (zz - mean) / std;
                (*n, zzz)
            })
            .collect();
        bar.finish();
        z
    }

    /// Генерирует ортотропное коррелированное узловое случайное поле с длиной корреляции
    /// в направлениях 2*lx, 2*ly, 2*lz
    pub fn generate_nodal_correlated_random_field_ortho(
        &self,
        lx: f32,
        ly: f32,
        lz: f32,
    ) -> HashMap<usize, f32> {
        let bar = ProgressBar::new(self.nodes.len() as u64).with_style(
            ProgressStyle::with_template(
                "{prefix} [{bar:40.cyan/blue}] {pos:>9}/{len:9} {percent:>3}% ETA {eta:>5}",
            )
            .unwrap()
            .progress_chars("=> "),
        );
        let l = lx.min(ly).min(lz);
        let rx = 3.0 * lx;
        let ry = 3.0 * ly;
        let rz = 3.0 * lz;
        let cell_size = 3.0 * l;
        let inv_2lx = 1.0 / 2.0 / lx / lx;
        let inv_2ly = 1.0 / 2.0 / ly / ly;
        let inv_2lz = 1.0 / 2.0 / lz / lz;
        let f: HashMap<usize, f32> = self
            .nodes
            .par_iter()
            .map(|(&n, _)| {
                let mut rng = rand::rng();
                (n, StandardNormal.sample(&mut rng))
            })
            .collect();
        let sh = self.nodes_spacial_hash(cell_size);
        let z: HashMap<usize, f32> = self
            .nodes
            .par_iter()
            .progress_with(bar.clone())
            .map(|(_, n0)| {
                let nds = self.node_neibors_ellipse(n0, &sh, cell_size, rx, ry, rz);
                let (s1, s2, s2q) =
                    nds.iter()
                        .fold((0.0f32, 0.0f32, 0.0f32), |(a1, a2, a3), &m| {
                            let n1 = &self.nodes[&m];
                            let dx2 = ((n1.x - n0.x) * (n1.x - n0.x)) as f32;
                            let dy2 = ((n1.y - n0.y) * (n1.y - n0.y)) as f32;
                            let dz2 = ((n1.z - n0.z) * (n1.z - n0.z)) as f32;
                            let w = fast_math::exp(-dx2 * inv_2lx - dy2 * inv_2ly - dz2 * inv_2lz);
                            (a1 + w * f[&m], a2 + w, a3 + w * w)
                        });
                let mean = s1 / s2;
                let sigma = (s2q / (s2 * s2)).sqrt();
                (n0.n, mean / sigma)
            })
            .collect();
        let zs = z.values().copied().collect::<Vec<f32>>();
        let (mean, std) = mean_and_std(&zs);
        let z: HashMap<usize, f32> = z
            .par_iter()
            .map(|(n, zz)| {
                let zzz = (zz - mean) / std;
                (*n, zzz)
            })
            .collect();
        bar.finish();
        z
    }

    /// Конвертирует узловые поля в ячеичные, осреднением по узлам элемента
    pub fn nodal_field_to_cell_field(&self, f: &HashMap<usize, f32>) -> HashMap<usize, f32> {
        self.solids
            .par_iter()
            .chain(self.shells.par_iter())
            .filter_map(|(&en, el)| {
                let n = el.nodes_count();
                if n == 0 {
                    return None;
                }
                let sum: f32 = (0..n).map(|nn| f[&el.nodes[nn]]).sum();
                Some((en, sum / n as f32))
            })
            .collect()
    }

    /// Ремапинг узлового поля на новую модель
    /// метод ближайшего соседа
    /// рекомендуется cell_size - размер наиболее грубой сетки из двух (исходной и на которую происходит перенос поля)
    pub fn remap_to_model_by_nearest(
        &self,
        cell_size: f32,
        other: &FEModel,
        nodal_field: &HashMap<usize, f32>,
    ) -> HashMap<usize, f32> {
        // let mut new_field: HashMap<usize, f32> = HashMap::with_capacity(other.nodes.len());
        let sh = self.nodes_spacial_hash(cell_size);
        other
            .nodes
            .par_iter()
            .map(|(_, v)| v)
            .filter_map(|n| {
                let nodes_d2 = self.node_neibors_with_distances2(n, &sh, cell_size, cell_size);
                let nn = nodes_d2
                    .iter()
                    .min_by(|a, b| a.1.total_cmp(b.1))
                    .map(|(&n, _)| n);
                if let Some(num) = nn {
                    if let Some(&val) = nodal_field.get(&num) {
                        Some((n.n, val))
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect()
    }

    /// Ремапинг узлового моля на новую модель
    /// с использованием метода обратного взвешивания расстояний (IDW)
    /// рекомендуется cell_size - размер наиболее грубой сетки из двух (исходной и на которую происходит перенос поля)
    pub fn remap_to_model_by_idw(
        &self,
        cell_size: f32,
        other: &FEModel,
        nodal_field: &HashMap<usize, f32>,
    ) -> HashMap<usize, f32> {
        // let mut new_field: HashMap<usize, f32> = HashMap::with_capacity(other.nodes.len());
        let sh = self.nodes_spacial_hash(cell_size);
        let tol2 = cell_size * cell_size / 1e6;
        other
            .nodes
            .par_iter()
            .map(|(_, v)| v)
            .map(|n| {
                let nodes_d2 = self.node_neibors_with_distances2(n, &sh, cell_size, cell_size);
                if let Some((&num, _)) = nodes_d2.iter().find(|(_, d2)| **d2 < tol2) {
                    return (n.n, nodal_field[&num]);
                }
                let (e, d) = nodes_d2.iter().fold((0.0, 0.0), |(e, d), (num, d2)| {
                    let w = 1.0 / d2 / d2;
                    (e + nodal_field[num] * w, d + w)
                });
                (n.n, e / d)
            })
            .collect()
    }

    pub fn generate_field_using_property_object(
        &self,
        fp: &FieldProperties,
    ) -> HashMap<usize, f32> {
        let normal_dist = match fp.field_type {
            FieldType::Isotropic(l) => self.generate_nodal_correlated_random_field(l / 2.0),
            FieldType::Orthotropic(lx, ly, lz) => {
                self.generate_nodal_correlated_random_field_ortho(lx / 2.0, ly / 2.0, lz / 2.0)
            }
        };
        norm_distribution_to_mean_std(&normal_dist, &fp.mean, &fp.std)
    }
}

/// Перевод нормального распределение со средним 0 и тандартным отклонением 1
/// в нормальное распределение со средним *mean* и абсолютным отклонением *std*
pub fn norm_distribution_to_mean_std(
    n: &HashMap<usize, f32>,
    mean: &f32,
    std: &f32,
) -> HashMap<usize, f32> {
    n.par_iter().map(|(&n, z)| (n, mean + std * z)).collect()
}

/// Стандартизация случайного поля к mean=0, std=1
pub fn standardize_field(f: &HashMap<usize, f32>) -> HashMap<usize, f32> {
    let v: Vec<f32> = f.values().copied().collect();
    let (mean, std) = mean_and_std(&v);
    if !(std > 0.0) {
        return f.clone();
    }
    f.par_iter().map(|(&k, &x)| (k, (x - mean) / std)).collect()
}

pub fn save_fields_to_abaqus(file_path: &str, fields: &Vec<&RandomField>) -> anyhow::Result<()> {
    let file = std::fs::File::create(file_path)?;
    let mut writer = std::io::BufWriter::new(file);
    writeln!(writer, "** Generated random fields")?;
    for &f in fields {
        let mut data: Vec<(usize, f32)> = f.data.clone().into_iter().collect();
        data.sort_by_key(|(k, _)| *k);

        writeln!(writer, "** {}", f)?;
        writeln!(
            writer,
            "*Initial conditions, type=FIELD, variable={}",
            &f.id
        )?;
        for (n, v) in data {
            writeln!(writer, "{}, {}", n, v)?;
        }
    }
    Ok(())
}

pub fn save_fields_to_lsdyna(file_path: &str, fields: &Vec<&RandomField>) -> anyhow::Result<()> {
    let n = fields.len();
    if n == 0 {
        return Err(anyhow::anyhow!("Пустой набор случайных полей"));
    }
    let file = std::fs::File::create(file_path)?;
    let mut writer = std::io::BufWriter::new(file);
    writeln!(writer, "*KEYWORD")?;
    for f in fields.iter() {
        writeln!(writer, "$ {}", f)?;
    }
    writeln!(writer, "*INITIAL_HISTORY_NODE")?;
    let mut node_numbers: Vec<usize> = fields[0].data.keys().cloned().collect();
    node_numbers.sort_unstable();
    for n_id in node_numbers {
        writeln!(writer, " {:>9} {:>9}", n_id, n)?;
        for i in 0..n {
            let rf = fields[i];
            if let Some(v) = rf.data.get(&n_id) {
                writeln!(writer, " {:>9} {:9.4e}", rf.id, v)?;
            }
        }
    }
    writeln!(writer, "*END")?;
    Ok(())
}

pub struct Perturbation {
    pub normal_field: HashMap<usize, f32>,
    pub field_type: FieldType,
}

impl Perturbation {
    pub fn new(model: &FEModel, field_type: FieldType) -> Self {
        match field_type {
            FieldType::Isotropic(lc) => Self {
                normal_field: model.generate_nodal_correlated_random_field(lc / 2.0),
                field_type: field_type,
            },
            FieldType::Orthotropic(lx, ly, lz) => Self {
                normal_field: model.generate_nodal_correlated_random_field_ortho(
                    lx / 2.0,
                    ly / 2.0,
                    lz / 2.0,
                ),
                field_type: field_type,
            },
        }
    }

    pub fn node_normalized_value(&self, nd_id: usize, mean: f32, std: f32) -> anyhow::Result<f32> {
        let v = self
            .normal_field
            .get(&nd_id)
            .context("В объекте perturbation нет узла с заданным идентификатором")?;
        Ok(mean + std * v)
    }

    pub fn node_ids(&self) -> Keys<'_, usize, f32> {
        self.normal_field.keys()
    }

    pub fn generate_for_norm_std(&self, mean: f32, std: f32) -> HashMap<usize, f32> {
        self.normal_field
            .iter()
            .map(|(&n, &v)| (n, mean + std * v))
            .collect()
    }
}

pub struct RandomField {
    pub id: usize,
    pub name: String,
    pub description: String,
    pub data: HashMap<usize, f32>,
}

impl RandomField {
    pub fn new() -> Self {
        Self {
            id: 1,
            name: "".to_string(),
            description: "".to_string(),
            data: HashMap::new(),
        }
    }

    pub fn with_id(mut self, id: usize) -> Self {
        self.id = id;
        self
    }

    pub fn with_name(mut self, name: String) -> Self {
        self.name = name;
        self
    }

    pub fn generate(mut self, p: &Perturbation, mean: f32, std: f32, inverse: bool) -> Self {
        let std = if inverse { -std } else { std };
        self.data = p.generate_for_norm_std(mean, std);
        self.description = format!(
            "mean: {}, std: {}, Field type: {:?}, inverse: {}",
            mean, std, p.field_type, inverse
        );
        self
    }
}

impl fmt::Display for RandomField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Random field {} ({}): {} values, {}",
            self.id,
            self.name,
            self.data.len(),
            self.description
        )
    }
}
