use crate::objects::femodel::FEModel;
use crate::utils::mean_and_std;
use fast_math;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rand_distr::{Distribution, StandardNormal};
use rayon::prelude::*;
use std::collections::HashMap;

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
}
