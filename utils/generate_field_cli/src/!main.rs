use clap::Parser;
use lsmesh::lsmodel::LsModel;
use std::path::PathBuf;

/// Генерация коррелированного случайного поля на КЭ-сетке.
#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    /// Путь к входному файлу сетки
    #[arg(short = 'i', long)]
    input: PathBuf,

    /// Масштабы корреляции: одно значение — изотропное поле,
    /// три значения (x,y,z) — ортотропное поле.
    /// Примеры: `-l 1.0` или `-s 1.0 2.0 0.5`
    #[arg(
        short = 'l',
        long,
        required = true,
        num_args = 1..=3,
        value_delimiter = ' ',
        value_name = "l",
        help = "Масштаб(ы) корреляции: 1 число (изотропия) или 3 числа (ортотропия)"
    )]
    correlation_scale: Vec<f32>,

    /// Среднее значение поля
    #[arg(short = 'm', long)]
    mean: f32,

    /// Стандартное отклонение поля
    #[arg(short = 'd', long)]
    std: f32,

    /// Путь к выходному VTK (по умолчанию: <input stem>.vtk)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() {
    let args = Args::parse();

    let output = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("vtk");
        p
    });

    let mut m = LsModel::default();
    println!("Чтение сетки");
    if let Some(file_name) = args.input.to_str() {
        m.read(file_name).expect("Ошибка чтения файла");
    }
    let l = args.correlation_scale;

    let z = match l.len() {
        1 => {
            println!("Генерация случайного нормального изотропного узлового поля");
            m.generate_nodal_correlated_random_field(l[0])
        }
        3 => {
            println!("Генерация случайного нормального ортотропного узлового поля");
            m.generate_nodal_correlated_random_field_ortho(l[0], l[1], l[2])
        }
        _ => {
            println!("Масштаб корреляции должен выражаться 1 или 3 числами");
            std::process::exit(2);
        }
    };

    println!("Пересчет в элементное нормальное поле");
    let cell_z = m.nodal_field_to_cell_field(&z);

    println!("Расчет значений поля с заданными параметрами распределения");
    let cell_val =
        lsmesh::lsmodel::utils::norm_distribution_to_mean_std(&cell_z, args.mean, args.std);

    println!("Сохранение файла vtk");
    if let Some(out_name) = output.to_str() {
        m.save_vtk(
            out_name,
            Some(&[("z", &z)]),
            Some(&[("ze", &cell_z), ("v", &cell_val)]),
        )
        .expect("Ошибка записи файла");
    }

    println!("Все операции завершены");
}
