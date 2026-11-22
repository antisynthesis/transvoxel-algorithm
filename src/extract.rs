//! Isosurface extraction: regular cells (modified Marching Cubes) and
//! transition cells, with the vertex-reuse scheme from Lengyel's dissertation.

use std::collections::HashMap;
use std::fmt;

use crate::block::{Block, Face, TransitionSides};
use crate::field::DensityField;
use crate::mesh::Mesh;
use crate::tables::{
    REGULAR_CELL_CLASS, REGULAR_CELL_DATA, REGULAR_VERTEX_DATA, TRANSITION_CASE_SAMPLES, TRANSITION_CELL_CLASS,
    TRANSITION_CELL_DATA, TRANSITION_CORNER_DATA, TRANSITION_VERTEX_DATA,
};

/// Parameters controlling extraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExtractOptions {
    /// Density value of the surface. Samples below it are inside.
    pub isovalue: f32,
    /// Thickness of the transition cells as a fraction of the block's cell size,
    /// in `[0, 1)`. Regular cells along transition faces are shrunk by this much.
    pub transition_width: f32,
}

impl Default for ExtractOptions {
    fn default() -> Self {
        Self { isovalue: 0.0, transition_width: 0.5 }
    }
}

impl ExtractOptions {
    /// Returns a copy with the given isovalue.
    pub fn with_isovalue(mut self, isovalue: f32) -> Self {
        self.isovalue = isovalue;
        self
    }

    /// Returns a copy with the given transition width.
    pub fn with_transition_width(mut self, transition_width: f32) -> Self {
        self.transition_width = transition_width;
        self
    }
}

/// Reasons extraction can be refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExtractError {
    /// The block has zero cells.
    EmptyBlock,
    /// The block's LOD is above 30, so its sample spacing would not fit in an `i32`.
    LodTooLarge,
    /// Transition faces were requested on an LOD-0 block. A transition face borders
    /// a block one level finer, and nothing is finer than LOD 0.
    TransitionAtLodZero,
    /// `transition_width` is not in `[0, 1)`.
    InvalidTransitionWidth,
    /// The block, plus its one-sample margin, does not fit in `i32` lattice coordinates.
    CoordinateOverflow,
}

impl fmt::Display for ExtractError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            ExtractError::EmptyBlock => "block has zero cells",
            ExtractError::LodTooLarge => "block level of detail is above 30",
            ExtractError::TransitionAtLodZero => "transition faces need a block with lod >= 1",
            ExtractError::InvalidTransitionWidth => "transition width must be in [0, 1)",
            ExtractError::CoordinateOverflow => "block coordinates overflow i32",
        })
    }
}

impl std::error::Error for ExtractError {}

/// Polygonizes `block` of `field` and returns the mesh.
///
/// Each face in `sides` must border a block whose LOD is exactly one lower
/// (twice the resolution); transition cells are generated there so that the two
/// meshes join without cracks.
pub fn extract<F: DensityField + ?Sized>(
    field: &F,
    block: &Block,
    sides: TransitionSides,
    options: &ExtractOptions,
) -> Result<Mesh, ExtractError> {
    let mut mesh = Mesh::new();
    extract_into(field, block, sides, options, &mut mesh)?;
    Ok(mesh)
}

/// Like [`extract`], but appends the geometry to an existing mesh.
///
/// Vertices are never shared with geometry already in `mesh`.
pub fn extract_into<F: DensityField + ?Sized>(
    field: &F,
    block: &Block,
    sides: TransitionSides,
    options: &ExtractOptions,
    mesh: &mut Mesh,
) -> Result<(), ExtractError> {
    validate(block, sides, options)?;
    let mut ex = Extraction::new(field, block, sides, options, mesh);
    ex.regular_cells();
    for face in sides.iter() {
        ex.transition_cells(face);
    }
    Ok(())
}

fn validate(block: &Block, sides: TransitionSides, options: &ExtractOptions) -> Result<(), ExtractError> {
    if block.cells == 0 {
        return Err(ExtractError::EmptyBlock);
    }
    if block.lod > 30 {
        return Err(ExtractError::LodTooLarge);
    }
    if !sides.is_empty() && block.lod == 0 {
        return Err(ExtractError::TransitionAtLodZero);
    }
    if !(0.0..1.0).contains(&options.transition_width) {
        return Err(ExtractError::InvalidTransitionWidth);
    }
    let step = 1i64 << block.lod;
    let extent = i64::from(block.cells) * step;
    for &o in &block.origin {
        let lo = i64::from(o) - step;
        let hi = i64::from(o) + extent + step;
        if lo < i64::from(i32::MIN) || hi > i64::from(i32::MAX) {
            return Err(ExtractError::CoordinateOverflow);
        }
    }
    Ok(())
}

// Under this crate's "below the isovalue is inside" convention the reference
// triangulations are already counter-clockwise seen from outside (for transition
// cells, with the face frames in `FACE_FRAME`), so no extra reversal is needed.
const REGULAR_REVERSE: bool = false;
const TRANSITION_REVERSE: bool = false;

const INVALID: u32 = u32::MAX;

/// For each face, the block axes along which the transition cell's `u` and `v`
/// sample axes run. Each frame satisfies `u × v = inward normal`, so all six
/// faces share one handedness and one winding rule.
const FACE_FRAME: [(usize, usize); 6] = [(1, 2), (2, 1), (2, 0), (0, 2), (0, 1), (1, 0)];

/// Where an output vertex lies, in some lattice's integer coordinates.
#[derive(Clone, Copy)]
enum At {
    /// Exactly on a lattice sample.
    Corner([i32; 3]),
    /// On the edge from the lower to the higher sample, at parameter `t` in `(0, 1)`.
    Edge([i32; 3], [i32; 3], f32),
}

impl At {
    fn new(a: [i32; 3], b: [i32; 3], t: f32) -> At {
        if t == 0.0 {
            At::Corner(a)
        } else if t == 1.0 {
            At::Corner(b)
        } else {
            At::Edge(a, b, t)
        }
    }
}

type EdgeKey = [i32; 6];

#[inline]
fn edge_key(a: [i32; 3], b: [i32; 3]) -> EdgeKey {
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    [a[0], a[1], a[2], b[0], b[1], b[2]]
}

/// Surface crossing parameter along an edge whose endpoints straddle the isovalue.
#[inline]
fn crossing(d0: f32, d1: f32, iso: f32) -> f32 {
    let t = (iso - d0) / (d1 - d0);
    if t.is_nan() {
        0.5
    } else {
        t.clamp(0.0, 1.0)
    }
}

#[inline]
fn to_f32(p: [i32; 3]) -> [f32; 3] {
    [p[0] as f32, p[1] as f32, p[2] as f32]
}

/// Interpolates between two lattice points. Every caller passes the endpoints in
/// the same (ascending) order, so neighbouring blocks produce bit-identical results.
#[inline]
fn lerp_point(a: [i32; 3], b: [i32; 3], t: f32) -> [f32; 3] {
    let (a, b) = (to_f32(a), to_f32(b));
    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])]
}

#[inline]
fn lerp3(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])]
}

#[inline]
fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 0.0 && len.is_finite() {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0; 3]
    }
}

#[inline]
fn add(a: [i32; 3], b: [i32; 3]) -> [i32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Offset of regular-cell corner `i`.
#[inline]
fn corner_offset(i: usize) -> [i32; 3] {
    [(i & 1) as i32, ((i >> 1) & 1) as i32, ((i >> 2) & 1) as i32]
}

/// `(u, v)` position of transition-cell sample `i` in half-cell units.
#[inline]
fn transition_sample_uv(i: usize) -> (i32, i32) {
    match i {
        0..=8 => ((i % 3) as i32, (i / 3) as i32),
        9 => (0, 0),
        10 => (2, 0),
        11 => (0, 2),
        _ => (2, 2),
    }
}

struct Extraction<'a, F: ?Sized> {
    field: &'a F,
    origin: [i32; 3],
    n: i32,
    step: i32,
    iso: f32,
    width: f32,
    sides: TransitionSides,
    /// Densities at block lattice points `-1..=n+1` on each axis.
    grid: Vec<f32>,
    gdim: usize,
    /// Regular vertices lying on a transition face, shared with the
    /// half-resolution face of the transition cells.
    shared_low: HashMap<EdgeKey, u32>,
    /// Full-resolution transition vertices on block edges, shared between the
    /// transition cells of two adjoining faces.
    shared_full: HashMap<EdgeKey, u32>,
    mesh: &'a mut Mesh,
}

impl<'a, F: DensityField + ?Sized> Extraction<'a, F> {
    fn new(field: &'a F, block: &Block, sides: TransitionSides, options: &ExtractOptions, mesh: &'a mut Mesh) -> Self {
        let n = block.cells as i32;
        let step = block.step();
        let origin = block.origin;
        let gdim = block.cells as usize + 3;
        let mut grid = Vec::with_capacity(gdim * gdim * gdim);
        for z in -1..=n + 1 {
            for y in -1..=n + 1 {
                for x in -1..=n + 1 {
                    grid.push(field.density(origin[0] + x * step, origin[1] + y * step, origin[2] + z * step));
                }
            }
        }
        Self {
            field,
            origin,
            n,
            step,
            iso: options.isovalue,
            width: options.transition_width,
            sides,
            grid,
            gdim,
            shared_low: HashMap::new(),
            shared_full: HashMap::new(),
            mesh,
        }
    }

    // ----------------------------------------------------------------------
    // Sampling
    // ----------------------------------------------------------------------

    /// Density at block lattice point `l` (each component in `-1..=n+1`).
    #[inline]
    fn sample(&self, l: [i32; 3]) -> f32 {
        let g = self.gdim;
        let i = |v: i32| (v + 1) as usize;
        self.grid[i(l[0]) + g * (i(l[1]) + g * i(l[2]))]
    }

    /// Central-difference gradient at block lattice point `l` (components in `0..=n`).
    #[inline]
    fn regular_gradient(&self, l: [i32; 3]) -> [f32; 3] {
        let mut g = [0.0; 3];
        for (a, ga) in g.iter_mut().enumerate() {
            let mut hi = l;
            let mut lo = l;
            hi[a] += 1;
            lo[a] -= 1;
            *ga = self.sample(hi) - self.sample(lo);
        }
        g
    }

    /// Central-difference gradient at LOD-0 point `p` with spacing `h`. Produces
    /// exactly what [`Self::regular_gradient`] produces in a block whose step is `h`.
    #[inline]
    fn field_gradient(&self, p: [i32; 3], h: i32) -> [f32; 3] {
        let mut g = [0.0; 3];
        for (a, ga) in g.iter_mut().enumerate() {
            let mut hi = p;
            let mut lo = p;
            hi[a] += h;
            lo[a] -= h;
            *ga = self.field.density(hi[0], hi[1], hi[2]) - self.field.density(lo[0], lo[1], lo[2]);
        }
        g
    }

    /// Block lattice point to LOD-0 lattice point.
    #[inline]
    fn world(&self, l: [i32; 3]) -> [i32; 3] {
        [self.origin[0] + l[0] * self.step, self.origin[1] + l[1] * self.step, self.origin[2] + l[2] * self.step]
    }

    // ----------------------------------------------------------------------
    // Vertex creation
    // ----------------------------------------------------------------------

    /// Creates a vertex of a regular cell (block lattice coordinates), applying the
    /// transition-face displacement where needed.
    fn make_regular_vertex(&mut self, at: At) -> u32 {
        let (la, lb, t) = match at {
            At::Corner(l) => (l, l, 0.0),
            At::Edge(a, b, t) => (a, b, t),
        };
        let (wa, wb) = (self.world(la), self.world(lb));
        let (mut position, normal) = match at {
            At::Corner(l) => (to_f32(wa), normalize(self.regular_gradient(l))),
            At::Edge(..) => {
                let g = lerp3(self.regular_gradient(la), self.regular_gradient(lb), t);
                (lerp_point(wa, wb, t), normalize(g))
            }
        };

        if !self.sides.is_empty() {
            if let Some(delta) = self.displacement(lerp_point(la, lb, t), normal) {
                let s = self.step as f32;
                for a in 0..3 {
                    position[a] += delta[a] * s;
                }
            }
        }

        let index = self.mesh.push_vertex(position, normal);
        if self.on_transition_face(la, lb) {
            self.shared_low.insert(edge_key(wa, wb), index);
        }
        index
    }

    /// The displacement (in cell units) that moves a vertex of a regular cell next
    /// to a transition face into the shrunken cell (dissertation section 4.4): the
    /// cell's boundary layer is compressed by `transition_width`, and the offset is
    /// projected onto the surface's tangent plane to preserve its shape.
    fn displacement(&self, local: [f32; 3], normal: [f32; 3]) -> Option<[f32; 3]> {
        let n = self.n as f32;
        let mut delta = [0.0f32; 3];
        let mut any = false;
        for a in 0..3 {
            if self.sides.contains(Face::ALL[2 * a]) && local[a] < 1.0 {
                delta[a] += (1.0 - local[a]) * self.width;
                any = true;
            }
            if self.sides.contains(Face::ALL[2 * a + 1]) && local[a] > n - 1.0 {
                delta[a] += (n - 1.0 - local[a]) * self.width;
                any = true;
            }
        }
        if !any {
            return None;
        }
        let nd = normal[0] * delta[0] + normal[1] * delta[1] + normal[2] * delta[2];
        Some([delta[0] - normal[0] * nd, delta[1] - normal[1] * nd, delta[2] - normal[2] * nd])
    }

    fn on_transition_face(&self, la: [i32; 3], lb: [i32; 3]) -> bool {
        self.sides.iter().any(|f| {
            let a = f.axis();
            let plane = if f.is_positive() { self.n } else { 0 };
            la[a] == plane && lb[a] == plane
        })
    }

    /// Creates a vertex on the full-resolution face of a transition cell
    /// (LOD-0 coordinates, sample spacing `h`). These vertices lie on the block
    /// boundary and are never displaced, so they match the finer neighbour exactly.
    fn make_full_vertex(&mut self, at: At, h: i32) -> u32 {
        let (position, normal) = match at {
            At::Corner(p) => (to_f32(p), normalize(self.field_gradient(p, h))),
            At::Edge(a, b, t) => {
                let g = lerp3(self.field_gradient(a, h), self.field_gradient(b, h), t);
                (lerp_point(a, b, t), normalize(g))
            }
        };
        self.mesh.push_vertex(position, normal)
    }

    /// `true` if the segment `a`–`b` (LOD-0 coordinates) lies on an edge of the
    /// block, where two faces meet.
    fn on_block_edge(&self, a: [i32; 3], b: [i32; 3]) -> bool {
        let e = self.n * self.step;
        let mut planes = 0;
        for k in 0..3 {
            if a[k] == b[k] && (a[k] == self.origin[k] || a[k] == self.origin[k] + e) {
                planes += 1;
            }
        }
        planes >= 2
    }

    // ----------------------------------------------------------------------
    // Regular cells
    // ----------------------------------------------------------------------

    /// Polygonizes every regular cell of the block.
    ///
    /// Vertex reuse follows section 3.3: each cell owns up to four vertices (the
    /// one on its maximal corner and those on the three edges ending there), and
    /// a cell finds a shared vertex in the preceding cell named by the vertex
    /// table's direction code. Two decks of per-cell slots (this z layer and the
    /// previous one) hold the vertices. The decks carry an extra row and column
    /// for the virtual cells just below the block's minimum faces, so vertices on
    /// those faces are shared too instead of duplicated.
    fn regular_cells(&mut self) {
        let n = self.n as usize;
        let dw = n + 1;
        let mut prev = vec![[INVALID; 4]; dw * dw];
        let mut cur = vec![[INVALID; 4]; dw * dw];
        let iso = self.iso;

        for z in 0..self.n {
            cur.fill([INVALID; 4]);
            for y in 0..self.n {
                for x in 0..self.n {
                    let cell = [x, y, z];
                    let mut density = [0.0f32; 8];
                    let mut case = 0usize;
                    for (i, d) in density.iter_mut().enumerate() {
                        *d = self.sample(add(cell, corner_offset(i)));
                        if *d < iso {
                            case |= 1 << i;
                        }
                    }
                    if case == 0 || case == 0xFF {
                        continue;
                    }

                    let data = &REGULAR_CELL_DATA[REGULAR_CELL_CLASS[case] as usize];
                    let mut local = [0u32; 12];
                    for (i, &code) in REGULAR_VERTEX_DATA[case][..data.vertex_count()].iter().enumerate() {
                        let v0 = ((code >> 4) & 0x0F) as usize;
                        let v1 = (code & 0x0F) as usize;
                        let t = crossing(density[v0], density[v1], iso);
                        let at = At::new(add(cell, corner_offset(v0)), add(cell, corner_offset(v1)), t);
                        let (dir, slot) = match at {
                            // A vertex on a corner is owned by the cell for which
                            // that corner is corner 7.
                            At::Corner(_) => {
                                let corner = if t == 0.0 { v0 } else { v1 };
                                ((corner ^ 7) as i32, 0)
                            }
                            At::Edge(..) => (((code >> 12) & 0x07) as i32, ((code >> 8) & 0x0F) as usize),
                        };
                        let ox = (x + 1 - (dir & 1)) as usize;
                        let oy = (y + 1 - ((dir >> 1) & 1)) as usize;
                        let deck = if dir & 4 != 0 { &mut prev } else { &mut cur };
                        let entry = &mut deck[ox + dw * oy][slot];
                        if *entry == INVALID {
                            *entry = self.make_regular_vertex(at);
                        }
                        local[i] = *entry;
                    }

                    for tri in data.indices().chunks_exact(3) {
                        let (a, b, c) = (local[tri[0] as usize], local[tri[1] as usize], local[tri[2] as usize]);
                        if REGULAR_REVERSE {
                            self.mesh.push_triangle(a, c, b);
                        } else {
                            self.mesh.push_triangle(a, b, c);
                        }
                    }
                }
            }
            std::mem::swap(&mut prev, &mut cur);
        }
    }

    // ----------------------------------------------------------------------
    // Transition cells
    // ----------------------------------------------------------------------

    /// Generates the transition cells along one face of the block (section 4.3).
    ///
    /// Each transition cell covers one boundary cell of the block. Its
    /// full-resolution face (samples `0..9`) lies on the block boundary and is
    /// sampled at half the block's spacing, matching the finer neighbour. Its
    /// half-resolution face (samples `9..=C`) coincides with the face of the
    /// shrunken regular cell, and its vertices are the regular cells' vertices.
    fn transition_cells(&mut self, face: Face) {
        let (ua, va) = FACE_FRAME[face as usize];
        let wa = face.axis();
        let n = self.n;
        let h = self.step / 2;
        let plane = if face.is_positive() { n } else { 0 };

        // LOD-0 position of the full-resolution sample at half-cell coordinates (hu, hv).
        let origin = self.origin;
        let half_point = move |hu: i32, hv: i32| -> [i32; 3] {
            let mut p = origin;
            p[ua] += hu * h;
            p[va] += hv * h;
            p[wa] += plane * 2 * h;
            p
        };
        // Block lattice position of the regular sample at cell coordinates (cu, cv).
        let lattice_point = move |cu: i32, cv: i32| -> [i32; 3] {
            let mut l = [0; 3];
            l[ua] = cu;
            l[va] = cv;
            l[wa] = plane;
            l
        };

        let sdim = (2 * n + 1) as usize;
        let mut samples = Vec::with_capacity(sdim * sdim);
        for hv in 0..=2 * n {
            for hu in 0..=2 * n {
                let p = half_point(hu, hv);
                samples.push(self.field.density(p[0], p[1], p[2]));
            }
        }

        let cw = (n + 1) as usize;
        let mut cache = vec![[INVALID; 10]; cw * cw];
        let iso = self.iso;

        for cv in 0..n {
            for cu in 0..n {
                let mut density = [0.0f32; 13];
                for (i, d) in density.iter_mut().enumerate().take(9) {
                    let (a, b) = transition_sample_uv(i);
                    *d = samples[(2 * cu + a) as usize + sdim * (2 * cv + b) as usize];
                }
                density[9] = density[0];
                density[10] = density[2];
                density[11] = density[6];
                density[12] = density[8];

                let mut case = 0usize;
                for (bit, &s) in TRANSITION_CASE_SAMPLES.iter().enumerate() {
                    if density[s] < iso {
                        case |= 1 << bit;
                    }
                }
                if case == 0 || case == 0x1FF {
                    continue;
                }

                let class = TRANSITION_CELL_CLASS[case];
                let data = &TRANSITION_CELL_DATA[(class & 0x7F) as usize];
                let inverted = class & 0x80 != 0;
                let mut local = [0u32; 12];
                // The center sample (4) is never shared with other cells, but
                // several edges of this cell may snap onto it.
                let mut center = INVALID;

                for (i, &code) in TRANSITION_VERTEX_DATA[case][..data.vertex_count()].iter().enumerate() {
                    let mut s0 = ((code >> 4) & 0x0F) as usize;
                    let mut s1 = (code & 0x0F) as usize;
                    let (mut uv0, mut uv1) = (transition_sample_uv(s0), transition_sample_uv(s1));
                    // Order the endpoints by ascending coordinate.
                    if uv0.0 + uv0.1 > uv1.0 + uv1.1 {
                        std::mem::swap(&mut s0, &mut s1);
                        std::mem::swap(&mut uv0, &mut uv1);
                    }
                    let t = crossing(density[s0], density[s1], iso);

                    local[i] = if s0 >= 9 {
                        // Half-resolution face: the vertex is the one the adjoining
                        // regular cell placed on the same edge.
                        let la = lattice_point(cu + uv0.0 / 2, cv + uv0.1 / 2);
                        let lb = lattice_point(cu + uv1.0 / 2, cv + uv1.1 / 2);
                        let at = At::new(la, lb, t);
                        let key = match at {
                            At::Corner(l) => edge_key(self.world(l), self.world(l)),
                            At::Edge(a, b, _) => edge_key(self.world(a), self.world(b)),
                        };
                        match self.shared_low.get(&key) {
                            Some(&index) => index,
                            None => self.make_regular_vertex(at),
                        }
                    } else {
                        let pa = half_point(2 * cu + uv0.0, 2 * cv + uv0.1);
                        let pb = half_point(2 * cu + uv1.0, 2 * cv + uv1.1);
                        let at = At::new(pa, pb, t);
                        let (dir, slot) = match at {
                            At::Corner(_) => {
                                let c = TRANSITION_CORNER_DATA[if t == 0.0 { s0 } else { s1 }];
                                ((c >> 4) as i32, (c & 0x0F) as usize)
                            }
                            At::Edge(..) => (((code >> 12) & 0x0F) as i32, ((code >> 8) & 0x0F) as usize),
                        };
                        if dir & 4 != 0 {
                            // Interior to this cell: never shared with other cells.
                            if matches!(at, At::Corner(_)) {
                                if center == INVALID {
                                    center = self.make_full_vertex(at, h);
                                }
                                center
                            } else {
                                self.make_full_vertex(at, h)
                            }
                        } else {
                            let ou = (cu + 1 - (dir & 1)) as usize;
                            let ov = (cv + 1 - ((dir >> 1) & 1)) as usize;
                            let entry = cache[ou + cw * ov][slot];
                            if entry != INVALID {
                                entry
                            } else {
                                let (ka, kb) = match at {
                                    At::Corner(p) => (p, p),
                                    At::Edge(a, b, _) => (a, b),
                                };
                                let index = if self.on_block_edge(ka, kb) {
                                    let key = edge_key(ka, kb);
                                    match self.shared_full.get(&key) {
                                        Some(&index) => index,
                                        None => {
                                            let index = self.make_full_vertex(at, h);
                                            self.shared_full.insert(key, index);
                                            index
                                        }
                                    }
                                } else {
                                    self.make_full_vertex(at, h)
                                };
                                cache[ou + cw * ov][slot] = index;
                                index
                            }
                        }
                    };
                }

                let reverse = inverted != TRANSITION_REVERSE;
                for tri in data.indices().chunks_exact(3) {
                    let (a, b, c) = (local[tri[0] as usize], local[tri[1] as usize], local[tri[2] as usize]);
                    if reverse {
                        self.mesh.push_triangle(a, c, b);
                    } else {
                        self.mesh.push_triangle(a, b, c);
                    }
                }
            }
        }
    }
}
