#[derive(Default, Debug)]
/// Object for rotation of vectors
pub struct Quarternion {
    data: [f64; 4],
    rotation_matrix: [[f64; 3]; 3],
    angle: f64,
    axis: Vector3D,
}

impl Quarternion {
    /// Create new quartenion for rotating vectors about *axis* on *angle*:
    ///
    /// **angle** in degrees;
    ///
    /// **axis** - axis of rotation Vector3D. It passes throw the origion (0,0,0)
    ///
    /// # Example
    ///
    /// ```
    /// use lsmodel::lsmesh::geom_utils::Quartenion;
    /// use lsmodel::lsmesh::geom_utils::Vector3D;
    /// use lsmodel::vec3d;
    /// let q = Quartenion(45, vec3d![1.0, 0.0, 0.0]);
    /// ```
    pub fn new(angle: f64, axis: Vector3D) -> Quarternion {
        let mut q = Quarternion {
            angle,
            axis: axis.normalized(),
            ..Quarternion::default()
        };
        let a = q.angle.to_radians();
        let a2 = a / 2.0;
        q.data[0] = a2.cos();
        let sina2 = a2.sin();
        let d = q.axis.coords();
        for i in 1..=3 {
            q.data[i] = d[i - 1] * sina2;
        }
        q.calculate_rotarion_matrix();
        q
    }

    fn calculate_rotarion_matrix(&mut self) {
        let w = self.data[0];
        let x = self.data[1];
        let y = self.data[2];
        let z = self.data[3];
        self.rotation_matrix[0][0] = 1.0 - 2.0 * y * y - 2.0 * z * z;
        self.rotation_matrix[0][1] = 2.0 * x * y - 2.0 * z * w;
        self.rotation_matrix[0][2] = 2.0 * x * z + 2.0 * y * w;
        self.rotation_matrix[1][0] = 2.0 * x * y + 2.0 * z * w;
        self.rotation_matrix[1][1] = 1.0 - 2.0 * x * x - 2.0 * z * z;
        self.rotation_matrix[1][2] = 2.0 * y * z - 2.0 * x * w;
        self.rotation_matrix[2][0] = 2.0 * x * z - 2.0 * y * w;
        self.rotation_matrix[2][1] = 2.0 * y * z + 2.0 * x * w;
        self.rotation_matrix[2][2] = 1.0 - 2.0 * x * x - 2.0 * y * y;
    }

    /// Rotates Vector3D *p* inplace
    ///
    /// # Example
    ///
    /// ```
    /// let q = Quarternion::new(180.0, vec3d![0.0, 10.0, 0.0]);
    /// let mut p = vec3d![15.0, 0.0, 0.0];
    /// q.rotate_point(&mut p);
    /// ```
    pub fn rotate_point(&self, p: &mut Vector3D) {
        let old_crds = p.coords();
        let mut_crds = p.coords_mut();
        for i in 0..3 {
            *mut_crds[i] = 0.0;
            for j in 0..3 {
                *mut_crds[i] += self.rotation_matrix[i][j] * old_crds[j];
            }
        }
    }

    /// Creates Vector3D by rotation of given Vector3D *p*
    ///
    /// # Example
    ///
    /// ```
    /// let q = Quarternion::new(180.0, vec3d![0.0, 10.0, 0.0]);
    /// let mut p = vec3d![15.0, 0.0, 0.0];
    /// let p2 = q.rotated_point(&p);
    /// ```
    pub fn rotated_point(&self, p: &Vector3D) -> Vector3D {
        let mut rp = Vector3D { ..*p };
        self.rotate_point(&mut rp);
        rp
    }
}

#[derive(Default, Debug)]
/// Object represented vector in 3D space defined by 3 spatial coordinates *x*, *y*, and *z*
///
/// # Example
///
/// ```
/// let v = Vector3D {x: 1.0, y: 0.0, z: 0.0};
/// let v = Vector3D::default();
/// let v = vec3d![1.0, 0.0, 0.0];
/// ```
pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3D {
    /// Returns the length of the Vector
    ///
    /// $\sqrt{x^2 + y^2 +z^2}$
    ///
    /// # Example
    ///
    /// ```
    /// let mut v = vec3d![2.0, 0.0, 0.0];
    /// let l: f64 = v.norm();
    /// assert_eq!(l, 1f64);
    /// assert_eq!(v.norm(), 1f64);
    /// ```
    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Normalizes vector implace. The length of vector after normalization is 1.0
    ///
    /// $ v = \frac{v}{|v|}$
    ///
    /// # Example
    ///
    /// ```
    /// let mut v = vec3d![4.0, 0.0, 0.0];
    /// v.normalize();
    /// assert_eq!(v.norm(), 1f64);
    /// ```
    pub fn normalize(&mut self) {
        let n = self.norm();
        self.x /= n;
        self.y /= n;
        self.z /= n;
    }

    /// Return normalized version of the vector
    pub fn normalized(&self) -> Vector3D {
        let mut p = Vector3D { ..*self };
        p.normalize();
        p
    }

    /// Scalar prouct of the vectors
    ///
    /// $ s = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z $
    ///
    /// # Example
    ///
    /// ```
    /// let v1 = vec3d![1.0, 0.0, 0.0];
    /// let v2 = vec3d![2.0, 1.0, 0.0];
    /// let s = v1.dot(&v2);
    /// ```
    pub fn dot(&self, other: &Vector3D) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Returns the cross product of two vectors.
    /// # Example
    ///
    /// ```
    /// let v1 = vec3d![1.0, 0.0, 0.0];
    /// let v2 = vec3d![2.0, 1.0, 0.0];
    /// let v: Vector3D = v1.cross(&v2);
    /// ```
    pub fn cross(&self, other: &Vector3D) -> Vector3D {
        let x = self.y * other.z - self.z * other.y;
        let y = self.z * other.x - self.x * other.z;
        let z = self.x * other.y - self.y * other.x;
        Vector3D { x, y, z }
    }

    /// Scales the vector inplace by scalar
    ///
    /// # Example
    ///
    /// ```
    /// let mut v1 = vec3d![1.0, 0.0, 0.0];
    /// v1.scale(10);
    /// assert_eq!(v1.norm(), 10f64);
    /// ```
    pub fn scale<T: Into<f64>>(&mut self, sc: T) {
        let sc: f64 = sc.into();
        self.x *= sc;
        self.y *= sc;
        self.z *= sc;
    }

    /// Returns the scaled version of the vector
    /// # Example
    ///
    /// ```
    /// let v1 = vec3d![1.0, 0.0, 0.0];
    /// let v2 = v1.scaled(10);
    /// assert_eq!(v2.norm(), 10f64);
    /// ```
    pub fn scaled<T: Into<f64>>(&self, sc: T) -> Vector3D {
        let mut p = Vector3D { ..*self };
        p.scale(sc);
        p
    }

    /// Returns shuufle product of three vectors
    ///
    /// v = v1.cross(v2).dot(v3)
    pub fn shuff_prod(&self, p2: &Vector3D, p3: &Vector3D) -> f64 {
        self.cross(p2).dot(p3)
    }

    /// Returns vector's coordinates as array [x, y, z]
    pub fn coords(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Returns vector's coordinates as mutable array [&mut x, &mut y, &mut z]
    pub fn coords_mut(&mut self) -> [&mut f64; 3] {
        [&mut self.x, &mut self.y, &mut self.z]
    }
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
        Vector3D {
            x: $x as f64,
            y: $y as f64,
            z: $z as f64,
        }
    };
}
