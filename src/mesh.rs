//! Triangle mesh output.

/// An indexed triangle mesh.
///
/// Positions are on the LOD-0 voxel lattice (one unit per voxel); use
/// [`Mesh::transform`] to move them into world space. Triangles are wound
/// counter-clockwise when viewed from outside the surface (right-handed
/// coordinates), and normals point outward.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    /// Vertex positions.
    pub positions: Vec<[f32; 3]>,
    /// Unit vertex normals, one per position.
    pub normals: Vec<[f32; 3]>,
    /// Triangle list, three indices per triangle.
    pub indices: Vec<u32>,
}

impl Mesh {
    /// Creates an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// `true` if the mesh has no triangles.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Removes all vertices and triangles, keeping the allocations.
    pub fn clear(&mut self) {
        self.positions.clear();
        self.normals.clear();
        self.indices.clear();
    }

    /// Iterates over triangles as index triples.
    pub fn triangles(&self) -> impl Iterator<Item = [u32; 3]> + '_ {
        self.indices.chunks_exact(3).map(|t| [t[0], t[1], t[2]])
    }

    /// Appends another mesh, offsetting its indices.
    pub fn append(&mut self, other: &Mesh) {
        let base = self.positions.len() as u32;
        self.positions.extend_from_slice(&other.positions);
        self.normals.extend_from_slice(&other.normals);
        self.indices.extend(other.indices.iter().map(|i| i + base));
    }

    /// Scales every position by `scale` and then adds `offset`.
    ///
    /// Normals are unaffected because the scale is uniform. A negative scale
    /// would mirror the mesh and is not supported.
    pub fn transform(&mut self, scale: f32, offset: [f32; 3]) {
        for p in &mut self.positions {
            for a in 0..3 {
                p[a] = p[a] * scale + offset[a];
            }
        }
    }

    pub(crate) fn push_vertex(&mut self, position: [f32; 3], normal: [f32; 3]) -> u32 {
        let index = self.positions.len() as u32;
        self.positions.push(position);
        self.normals.push(normal);
        index
    }

    #[inline]
    pub(crate) fn push_triangle(&mut self, a: u32, b: u32, c: u32) {
        // Degenerate triangles appear when several vertices snap to the same
        // voxel corner (the density equals the isovalue exactly); drop them.
        if a != b && b != c && c != a {
            self.indices.extend_from_slice(&[a, b, c]);
        }
    }
}
