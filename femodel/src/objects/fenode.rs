use crate::utils::HasNum;
use crate::utils::parse_sub_str;
use std::fmt::Display;
use std::fmt::Write;

#[derive(Debug, Default)]
/// Объект узла конечно-элементной сетки
#[derive(Clone)]
pub struct FENode {
    /// Номер узла
    pub n: usize,
    /// Координата X
    pub x: f64,
    /// Координата Y
    pub y: f64,
    /// Координата Z
    pub z: f64,
}

impl Display for FENode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Node {}: ({}, {}, {})", self.n, self.x, self.y, self.z)
    }
}

impl FENode {
    /// Перемещение узла на вектор u
    pub fn translate(&mut self, u: &[f64]) {
        self.x += u.first().unwrap_or(&0.0);
        self.y += u.get(1).unwrap_or(&0.0);
        self.z += u.get(2).unwrap_or(&0.0);
    }
    /// Координаты узла как неизменяемый массив (только для чтения)
    pub fn coordinates(&self) -> [&f64; 3] {
        [&self.x, &self.y, &self.z]
    }
    /// Координаты узла как массив с возможностью изменения значений
    pub fn coordinates_mut(&mut self) -> [&mut f64; 3] {
        [&mut self.x, &mut self.y, &mut self.z]
    }
    /// Расстояние до другого узла
    pub fn distance_to_node(&self, nd2: &FENode) -> f64 {
        f64::sqrt((self.x - nd2.x).powi(2) + (self.y - nd2.y).powi(2) + (self.z - nd2.z).powi(2))
    }

    /// Квадрат расстояния до другого узла
    pub fn distance_squared_to_node(&self, nd2: &FENode) -> f32 {
        ((self.x - nd2.x).powi(2) + (self.y - nd2.y).powi(2) + (self.z - nd2.z).powi(2)) as f32
    }

    /// возвращает индекс ячейки размером cell_size, в которую попадает узел в пространстве
    pub fn spacial_hash(&self, cell_size: f32) -> (isize, isize, isize) {
        (
            (self.x / cell_size as f64).floor() as isize,
            (self.y / cell_size as f64).floor() as isize,
            (self.z / cell_size as f64).floor() as isize,
        )
    }
}

impl From<&str> for FENode {
    /// Инициализация узла из строки.
    /// Формат строки: N, X, Y, Z
    /// Числа могут разделяться запятыми или используется формат карты *NODE:
    /// 8 позиций под номер и по 16 позиций под координаты
    fn from(value: &str) -> Self {
        // Если используется разделение запятыми, по парсинг данных происходит здесь
        if value.as_bytes().contains(&b',') {
            let v: Vec<_> = value.split(',').map(|substr| substr.trim()).collect();
            let n = v
                .first()
                .unwrap_or(&"")
                .parse::<usize>()
                .unwrap_or_default();
            let x = v.get(1).unwrap_or(&"").parse::<f64>().unwrap_or_default();
            let y = v.get(2).unwrap_or(&"").parse::<f64>().unwrap_or_default();
            let z = v.get(3).unwrap_or(&"").parse::<f64>().unwrap_or_default();
            FENode { n, x, y, z }
        // Иначе данные парсятся здесь
        } else {
            let n: usize = parse_sub_str(value, 0, 8);
            let x: f64 = parse_sub_str(value, 8, 24);
            let y: f64 = parse_sub_str(value, 24, 40);
            let z: f64 = parse_sub_str(value, 40, 56);
            FENode { n, x, y, z }
        }
    }
}

impl From<&FENode> for String {
    fn from(nd: &FENode) -> String {
        let mut buf: String = String::with_capacity(56);
        write!(
            &mut buf,
            " {:7} {:15.6e} {:15.6e} {:15.6e}",
            nd.n, nd.x, nd.y, nd.z
        )
        .unwrap();
        buf
    }
}

impl HasNum for FENode {
    fn mut_num(&mut self) -> &mut usize {
        &mut self.n
    }
}
