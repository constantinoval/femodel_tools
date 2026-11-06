pub mod objects;
pub mod readers;
pub mod utils;

pub mod prelude {
    pub use crate::objects::feelement::FEElement;
    pub use crate::objects::femodel::FEModel;
    pub use crate::objects::femodel::{CellScalarField, NodeScalarField};
    pub use crate::objects::fenode::FENode;
}

#[macro_export]
/// Macro for creation of Vector3D
///
/// # Example
///
/// ```
/// let v = vec3d![1.0, 0.0, 0.0];
/// let v = Vector3D {x: 1.0, y: 0.0, z: 0.0};
/// ```
macro_rules! vec3d {
    [$x: expr, $y: expr, $z: expr] => {
        femodel::utils::geom_utils::Vector3D {
            x: $x as f64,
            y: $y as f64,
            z: $z as f64,
        }
    };
}
