//! Blocks, faces and transition-side sets.

/// One of the six faces of a [`Block`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Face {
    /// The face at the minimum x coordinate.
    NegX = 0,
    /// The face at the maximum x coordinate.
    PosX = 1,
    /// The face at the minimum y coordinate.
    NegY = 2,
    /// The face at the maximum y coordinate.
    PosY = 3,
    /// The face at the minimum z coordinate.
    NegZ = 4,
    /// The face at the maximum z coordinate.
    PosZ = 5,
}

impl Face {
    /// All six faces, in declaration order.
    pub const ALL: [Face; 6] = [Face::NegX, Face::PosX, Face::NegY, Face::PosY, Face::NegZ, Face::PosZ];

    /// The axis perpendicular to this face (`0` = x, `1` = y, `2` = z).
    #[inline]
    pub const fn axis(self) -> usize {
        self as usize / 2
    }

    /// `true` for the faces at the maximum coordinate of their axis.
    #[inline]
    pub const fn is_positive(self) -> bool {
        self as usize % 2 == 1
    }

    /// The face on the opposite side of the block.
    #[inline]
    pub const fn opposite(self) -> Face {
        Face::ALL[self as usize ^ 1]
    }

    #[inline]
    const fn bit(self) -> u8 {
        1 << self as u8
    }
}

/// A set of block faces that border a block of *higher* resolution (lower LOD).
///
/// Transition cells are generated along every face in the set, and the regular
/// cells next to those faces are shrunk to make room for them.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct TransitionSides(u8);

impl TransitionSides {
    /// No transition faces.
    pub const NONE: Self = Self(0);
    /// All six faces.
    pub const ALL: Self = Self(0x3F);

    /// Creates an empty set.
    #[inline]
    pub const fn empty() -> Self {
        Self::NONE
    }

    /// Returns `true` if `face` is in the set.
    #[inline]
    pub const fn contains(self, face: Face) -> bool {
        self.0 & face.bit() != 0
    }

    /// Returns `true` if the set is empty.
    #[inline]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Returns a copy of the set with `face` added.
    #[inline]
    pub const fn with(self, face: Face) -> Self {
        Self(self.0 | face.bit())
    }

    /// Returns a copy of the set with `face` removed.
    #[inline]
    pub const fn without(self, face: Face) -> Self {
        Self(self.0 & !face.bit())
    }

    /// Adds `face` to the set.
    #[inline]
    pub fn insert(&mut self, face: Face) {
        self.0 |= face.bit();
    }

    /// Removes `face` from the set.
    #[inline]
    pub fn remove(&mut self, face: Face) {
        self.0 &= !face.bit();
    }

    /// Iterates over the faces in the set.
    pub fn iter(self) -> impl Iterator<Item = Face> {
        Face::ALL.into_iter().filter(move |f| self.contains(*f))
    }

    /// The raw bit mask (bit `i` = `Face::ALL[i]`).
    #[inline]
    pub const fn bits(self) -> u8 {
        self.0
    }

    /// Builds a set from a raw bit mask; bits above the sixth are ignored.
    #[inline]
    pub const fn from_bits_truncate(bits: u8) -> Self {
        Self(bits & 0x3F)
    }
}

impl From<Face> for TransitionSides {
    fn from(face: Face) -> Self {
        Self(face.bit())
    }
}

impl<const N: usize> From<[Face; N]> for TransitionSides {
    fn from(faces: [Face; N]) -> Self {
        faces.into_iter().collect()
    }
}

impl FromIterator<Face> for TransitionSides {
    fn from_iter<I: IntoIterator<Item = Face>>(iter: I) -> Self {
        iter.into_iter().fold(Self::NONE, Self::with)
    }
}

impl core::ops::BitOr for TransitionSides {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl core::ops::BitOr<Face> for TransitionSides {
    type Output = Self;
    fn bitor(self, rhs: Face) -> Self {
        self.with(rhs)
    }
}

impl core::ops::BitOr for Face {
    type Output = TransitionSides;
    fn bitor(self, rhs: Face) -> TransitionSides {
        TransitionSides::from(self).with(rhs)
    }
}

impl core::ops::BitOrAssign for TransitionSides {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl core::ops::BitOrAssign<Face> for TransitionSides {
    fn bitor_assign(&mut self, rhs: Face) {
        self.insert(rhs);
    }
}

/// A cubic block of cells to polygonize.
///
/// All coordinates are expressed on the *finest* voxel lattice (LOD 0). A block
/// at level of detail `lod` samples every `2^lod` voxels, so it spans
/// `cells * 2^lod` voxels along each axis, starting at `origin`.
///
/// Two adjacent blocks fit together without cracks when they have the same LOD,
/// or when the coarser block (higher `lod`, exactly one level higher) lists the
/// shared face in its [`TransitionSides`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Block {
    /// Minimum corner of the block on the LOD-0 lattice.
    pub origin: [i32; 3],
    /// Number of cells along each axis.
    pub cells: u32,
    /// Level of detail: the block samples every `2^lod` voxels.
    pub lod: u8,
}

impl Block {
    /// Creates a block.
    #[inline]
    pub const fn new(origin: [i32; 3], cells: u32, lod: u8) -> Self {
        Self { origin, cells, lod }
    }

    /// Distance between neighbouring samples on the LOD-0 lattice (`2^lod`).
    #[inline]
    pub const fn step(&self) -> i32 {
        1 << self.lod
    }

    /// Edge length of the block on the LOD-0 lattice (`cells * 2^lod`).
    #[inline]
    pub const fn extent(&self) -> i32 {
        self.cells as i32 * self.step()
    }

    /// Maximum corner of the block on the LOD-0 lattice.
    #[inline]
    pub const fn max(&self) -> [i32; 3] {
        let e = self.extent();
        [self.origin[0] + e, self.origin[1] + e, self.origin[2] + e]
    }
}
