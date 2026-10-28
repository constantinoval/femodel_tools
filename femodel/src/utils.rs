use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::str::FromStr;

pub fn parse_sub_str<T>(s: &str, start: u8, end: u8) -> T
where
    T: FromStr + Default,
{
    s.get(start as usize..end as usize)
        .unwrap_or("")
        .trim()
        .parse::<T>()
        .unwrap_or_default()
}

pub type Matrix3x3 = [[f64; 3]; 3];

pub fn determinant(m: &Matrix3x3) -> f64 {
    let mut result: f64 = 0.0;
    result += m[0][0] * m[1][1] * m[2][2];
    result += m[0][1] * m[1][2] * m[2][0];
    result += m[1][0] * m[2][1] * m[0][2];
    result -= m[0][2] * m[1][1] * m[2][0];
    result -= m[0][0] * m[1][2] * m[2][1];
    result -= m[2][2] * m[1][0] * m[0][1];
    result
}

pub trait HasNum {
    fn mut_num(&mut self) -> &mut usize;
}

pub fn renumber_map<T>(
    m: BTreeMap<usize, T>,
    start_number: usize,
    return_new_keys: bool,
) -> (BTreeMap<usize, T>, HashMap<usize, usize>)
where
    T: HasNum,
{
    let new_keys: HashMap<usize, usize> = if return_new_keys {
        m.keys()
            .enumerate()
            .map(|(i, &k)| (k, i + start_number))
            .collect()
    } else {
        HashMap::new()
    };
    let mut m: BTreeMap<usize, T> = m
        .into_values()
        .enumerate()
        .map(|(i, v)| (i + start_number, v))
        .collect();
    m.iter_mut().for_each(|(k, v)| *v.mut_num() = *k);
    (m, new_keys)
}

/// Вычисляет среднее значение и выборочное стандартное отклонение (std) для среза f64.
pub fn mean_and_std(data: &[f32]) -> (f32, f32) {
    let n = data.len() as f32;
    if n < 2.0 {
        return (0.0, 0.0);
    }
    let mean = data.par_iter().map(|&x| x).sum::<f32>() / n;
    let var = data
        .par_iter()
        .map(|&x| (x - mean) * (x - mean))
        .sum::<f32>()
        / (n - 1.0);
    (mean, var.max(0.0).sqrt())
}

pub fn spacial_hash(crds: &[f64; 3], cell_size: f64) -> (isize, isize, isize) {
    (
        (crds[0] / cell_size).floor() as isize,
        (crds[1] / cell_size).floor() as isize,
        (crds[2] / cell_size).floor() as isize,
    )
}

/// Перевод нормального распределение со средним 0 и тандартным отклонением 1
/// в нормальное распределение со средним *mean* и абсолютным отклонением *std*
pub fn norm_distribution_to_mean_std(
    n: &HashMap<usize, f32>,
    mean: f32,
    std: f32,
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

pub mod generations;
pub mod geom_utils;
