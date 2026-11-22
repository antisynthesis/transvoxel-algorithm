# transvoxel-algorithm

A complete Rust implementation of Eric Lengyel's [Transvoxel Algorithm](https://transvoxel.org/):
it turns voxel data into triangle meshes at multiple levels of detail (LOD) with no
cracks or seams between blocks of different resolution.

- **Regular cells**: modified Marching Cubes that resolves ambiguous faces consistently.
- **Transition cells**: generated on any block face that borders a block of twice the resolution.
- **Vertex reuse**: the reference scheme from the dissertation, so each shared vertex is created
  once and indexed. Vertices on the block's minimum faces are shared too.
- Vertices placed exactly on voxel corners when a density equals the isovalue.
- Gradient (central-difference) normals.
- Regular vertices next to transition faces are displaced, with the offset projected onto the
  tangent plane, as in section 4.4 of the dissertation.
- The official lookup tables are exposed in `transvoxel::tables`.
- No dependencies and no `unsafe`.

## Usage

```toml
[dependencies]
transvoxel-algorithm = { git = "https://github.com/antisynthesis/transvoxel-algorithm" }
```

The library is imported as `transvoxel`:

```rust
use transvoxel::{extract, Block, ExtractOptions, Face, TransitionSides};

// Any `Fn(i32, i32, i32) -> f32` is a density field. Negative = inside.
let sphere = |x: i32, y: i32, z: i32| {
    let (x, y, z) = (x as f32 - 32.0, y as f32 - 32.0, z as f32 - 32.0);
    (x * x + y * y + z * z).sqrt() - 20.0
};

// 16³ cells at LOD 1 (samples every 2nd voxel), covering voxels 0..32.
let block = Block::new([0, 0, 0], 16, 1);

// The +X neighbour is at LOD 0, so stitch that face with transition cells.
let mesh = extract(&sphere, &block, Face::PosX.into(), &ExtractOptions::default()).unwrap();

for [a, b, c] in mesh.triangles() {
    let _ = (mesh.positions[a as usize], mesh.normals[a as usize], b, c);
}
```

### Concepts

| Item | Meaning |
| --- | --- |
| `DensityField` | Scalar field on the integer LOD-0 voxel lattice. Implemented for closures and for `VoxelGrid`. |
| `Block { origin, cells, lod }` | A cube of `cells³` cells that samples every `2^lod` voxels, starting at `origin` (LOD-0 coordinates). |
| `TransitionSides` | The faces of a block that border a block exactly **one LOD finer**. |
| `ExtractOptions` | `isovalue` (default `0.0`) and `transition_width` (default `0.5`, as a fraction of a cell). |
| `Mesh` | `positions`, `normals` and `indices` (`u32` triangle list). Use `extract_into` to append several blocks into one mesh. |

### Conventions

- Densities **below** the isovalue are inside (solid), as in signed distance fields.
- Positions are in LOD-0 voxel units. `Mesh::transform(scale, offset)` maps them to world space.
- Triangles are counter-clockwise seen from outside, and normals point outward (right-handed).
- Blocks need the field to be defined one sample beyond each face, for normals.

### Choosing transition sides

For each face of a block, look at the neighbouring block across it:

- same LOD: nothing to do;
- one LOD **finer** (half the step): add that face to the block's `TransitionSides`;
- one LOD **coarser**: nothing to do here. The coarser block adds the transition face on its side.

Neighbouring LODs must differ by at most one level, as with any restricted octree.

## Example

```sh
cargo run --release --example terrain -- terrain.obj
```

builds a two-level terrain with a cave, made of 26 LOD-1 blocks around 8 LOD-0 blocks, and
writes it as an OBJ file.

## Testing

`cargo test` checks the structural invariants of the tables, and builds scenes in which a coarse
block is surrounded on all six faces by finer blocks. The surface crosses the coarse block's
faces, edges and corners at LOD pairs 1/0 and 2/1, with several transition widths. The tests
then verify that the combined mesh is closed and consistently oriented: every directed edge is
matched by exactly one reverse edge.

## Known limits

These come from the algorithm itself:

- With small `transition_width` values on coarse blocks, the tangent-plane projection can slide a
  vertex further than the transition slab is deep. A few thin slivers then fold over. The mesh
  stays closed.
- Where a transition face meets a face bordering a block of the **same** LOD, the projected
  displacement can move boundary vertices slightly off that shared face.

## License

MIT. The lookup tables in `src/tables.rs` are © 2009 Eric Lengyel, MIT licensed
(<https://github.com/EricLengyel/Transvoxel>), and the notice is kept in that file.
