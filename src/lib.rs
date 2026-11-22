//! # Transvoxel
//!
//! A complete Rust implementation of Eric Lengyel's
//! [Transvoxel Algorithm](https://transvoxel.org/) for extracting triangle meshes
//! from voxel data at multiple levels of detail without cracks or seams.
//!
//! The algorithm has two parts:
//!
//! * **Regular cells** use a modified Marching Cubes whose 18 equivalence
//!   classes resolve ambiguous faces consistently, so neighbouring cells and
//!   blocks always agree on the surface.
//! * **Transition cells** sit on a block face that borders a block of twice the
//!   resolution. Their full-resolution side matches the finer block exactly, and
//!   their half-resolution side matches the coarser block's regular cells, which
//!   are shrunk to make room for them.
//!
//! This crate implements both, including the reference vertex-reuse scheme
//! (shared vertices are created once and indexed), vertex placement exactly on
//! voxel corners, interpolated gradient normals, and the projected displacement
//! of regular vertices next to transition faces.
//!
//! ## Example
//!
//! ```
//! use transvoxel::{extract, Block, ExtractOptions, Face, TransitionSides};
//!
//! // A sphere of radius 20 voxels centred at (32, 32, 32); negative inside.
//! let sphere = |x: i32, y: i32, z: i32| {
//!     let (x, y, z) = (x as f32 - 32.0, y as f32 - 32.0, z as f32 - 32.0);
//!     (x * x + y * y + z * z).sqrt() - 20.0
//! };
//!
//! // A 16³-cell block at LOD 1 (every 2nd voxel) covering voxels 0..32 on each axis.
//! let block = Block::new([0, 0, 0], 16, 1);
//!
//! // Its +X face borders a block at LOD 0, so stitch that face with transition cells.
//! let sides = TransitionSides::from(Face::PosX);
//!
//! let mesh = extract(&sphere, &block, sides, &ExtractOptions::default()).unwrap();
//! assert!(mesh.triangle_count() > 0);
//! assert_eq!(mesh.positions.len(), mesh.normals.len());
//! ```
//!
//! ## Conventions
//!
//! * Densities **below** [`ExtractOptions::isovalue`] are inside the surface, as
//!   in signed distance fields.
//! * Positions are on the LOD-0 voxel lattice; see [`Mesh::transform`].
//! * Triangles are counter-clockwise seen from outside, and normals point outward.
//! * A face listed in [`TransitionSides`] must border a block exactly one LOD finer.
//!   Faces bordering blocks of the same LOD need nothing special.
//!
//! The raw lookup tables are available in [`tables`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod block;
mod extract;
mod field;
mod mesh;
pub mod tables;

pub use block::{Block, Face, TransitionSides};
pub use extract::{extract, extract_into, ExtractError, ExtractOptions};
pub use field::{DensityField, VoxelGrid};
pub use mesh::Mesh;
