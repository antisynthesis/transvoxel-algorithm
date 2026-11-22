//! Voxel data sources.

/// A scalar field sampled on the integer LOD-0 voxel lattice.
///
/// Points whose density is **below** the isovalue are *inside* (solid) and points
/// at or above it are *outside* (empty). This matches signed distance fields,
/// which are negative inside the shape. Normals point from inside to outside.
///
/// The extractor samples a margin around each block (one sample beyond each
/// face, for gradients), so the field must be defined outside the block too.
pub trait DensityField {
    /// The density at lattice point `(x, y, z)`.
    fn density(&self, x: i32, y: i32, z: i32) -> f32;
}

impl<F> DensityField for F
where
    F: Fn(i32, i32, i32) -> f32,
{
    #[inline]
    fn density(&self, x: i32, y: i32, z: i32) -> f32 {
        self(x, y, z)
    }
}

/// A dense, axis-aligned box of voxel samples.
///
/// Lookups outside the box are clamped to the nearest sample on its border.
#[derive(Clone, Debug, PartialEq)]
pub struct VoxelGrid {
    origin: [i32; 3],
    dims: [usize; 3],
    data: Vec<f32>,
}

impl VoxelGrid {
    /// Creates a grid of `dims` samples starting at `origin`, filled with `value`.
    ///
    /// # Panics
    /// Panics if any dimension is zero.
    pub fn new(origin: [i32; 3], dims: [usize; 3], value: f32) -> Self {
        assert!(dims.iter().all(|&d| d > 0), "VoxelGrid dimensions must be non-zero");
        Self { origin, dims, data: vec![value; dims[0] * dims[1] * dims[2]] }
    }

    /// Creates a grid by evaluating `f` at every lattice point of the box.
    pub fn from_fn(origin: [i32; 3], dims: [usize; 3], mut f: impl FnMut(i32, i32, i32) -> f32) -> Self {
        let mut grid = Self::new(origin, dims, 0.0);
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let i = grid.linear(x, y, z);
                    grid.data[i] = f(origin[0] + x as i32, origin[1] + y as i32, origin[2] + z as i32);
                }
            }
        }
        grid
    }

    /// The minimum lattice point covered by the grid.
    pub fn origin(&self) -> [i32; 3] {
        self.origin
    }

    /// The number of samples along each axis.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// The raw samples, x-major (x varies fastest, then y, then z).
    pub fn data(&self) -> &[f32] {
        &self.data
    }

    /// Mutable access to the raw samples.
    pub fn data_mut(&mut self) -> &mut [f32] {
        &mut self.data
    }

    #[inline]
    fn linear(&self, x: usize, y: usize, z: usize) -> usize {
        x + self.dims[0] * (y + self.dims[1] * z)
    }

    #[inline]
    fn local(&self, p: [i32; 3]) -> Option<[usize; 3]> {
        let mut out = [0; 3];
        for a in 0..3 {
            let v = p[a].checked_sub(self.origin[a])?;
            if v < 0 || v as usize >= self.dims[a] {
                return None;
            }
            out[a] = v as usize;
        }
        Some(out)
    }

    /// The sample at `p`, or `None` if `p` lies outside the grid.
    pub fn get(&self, p: [i32; 3]) -> Option<f32> {
        self.local(p).map(|[x, y, z]| self.data[self.linear(x, y, z)])
    }

    /// Sets the sample at `p`. Returns `false` (and does nothing) if `p` lies outside the grid.
    pub fn set(&mut self, p: [i32; 3], value: f32) -> bool {
        match self.local(p) {
            Some([x, y, z]) => {
                let i = self.linear(x, y, z);
                self.data[i] = value;
                true
            }
            None => false,
        }
    }
}

impl DensityField for VoxelGrid {
    #[inline]
    fn density(&self, x: i32, y: i32, z: i32) -> f32 {
        let c = |v: i32, a: usize| -> usize {
            let l = i64::from(v) - i64::from(self.origin[a]);
            l.clamp(0, self.dims[a] as i64 - 1) as usize
        };
        self.data[self.linear(c(x, 0), c(y, 1), c(z, 2))]
    }
}
