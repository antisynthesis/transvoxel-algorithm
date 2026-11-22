//! End-to-end tests: closed, consistently oriented meshes across blocks of
//! different levels of detail.

mod common;

use common::*;
use std::collections::HashSet;
use transvoxel::*;

/// Every block (of `cells` cells at `lod`) tiling the box `[lo, hi)³`, except
/// those inside `hole`.
fn tiles(lo: i32, hi: i32, cells: u32, lod: u8, hole: Option<([i32; 3], [i32; 3])>) -> Vec<Block> {
    let e = cells as i32 * (1 << lod);
    let mut out = Vec::new();
    for z in (lo..hi).step_by(e as usize) {
        for y in (lo..hi).step_by(e as usize) {
            for x in (lo..hi).step_by(e as usize) {
                let o = [x, y, z];
                if let Some((a, b)) = hole {
                    if (0..3).all(|k| o[k] >= a[k] && o[k] < b[k]) {
                        continue;
                    }
                }
                out.push(Block::new(o, cells, lod));
            }
        }
    }
    out
}

fn build(field: &impl DensityField, blocks: &[(Block, TransitionSides)], options: &ExtractOptions) -> Mesh {
    let mut mesh = Mesh::new();
    for (block, sides) in blocks {
        extract_into(field, block, *sides, options, &mut mesh).unwrap();
    }
    mesh
}

fn assert_closed_outward(mesh: &Mesh, center: [f32; 3]) {
    assert!(mesh.triangle_count() > 0);
    if let Err(e) = check_closed(mesh) {
        panic!("mesh is not closed: {e}");
    }
    assert!(signed_volume(mesh) > 0.0, "mesh is inside out");
    let f = outward_fraction(mesh, center);
    assert!(f > 0.97, "only {f} of triangles face outward");
}

#[test]
fn single_block_sphere_is_closed() {
    let c = [16.0, 16.0, 16.0];
    let field = sphere(c, 10.3);
    let mesh =
        extract(&field, &Block::new([0, 0, 0], 32, 0), TransitionSides::NONE, &ExtractOptions::default()).unwrap();
    assert_closed_outward(&mesh, c);

    // Volume close to the analytic value.
    let v = signed_volume(&mesh);
    let exact = 4.0 / 3.0 * std::f64::consts::PI * 10.3f64.powi(3);
    assert!((v - exact).abs() / exact < 0.02, "volume {v} vs {exact}");

    // Vertex reuse: no two vertices share a position, and the surface is a sphere (χ = 2).
    let unique: HashSet<_> = mesh.positions.iter().map(|p| p.map(f32::to_bits)).collect();
    assert_eq!(unique.len(), mesh.vertex_count());
    let mut edges = HashSet::new();
    for t in mesh.triangles() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            edges.insert((a.min(b), a.max(b)));
        }
    }
    let chi = mesh.vertex_count() as i64 - edges.len() as i64 + mesh.triangle_count() as i64;
    assert_eq!(chi, 2);

    // Normals point outward and are unit length.
    for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
        let r = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
        let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        let dot = (r[0] * n[0] + r[1] * n[1] + r[2] * n[2]) / len;
        assert!(dot > 0.95, "normal {n:?} at {p:?}");
        assert!(((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) - 1.0).abs() < 1e-4);
    }
}

#[test]
fn same_lod_blocks_join_without_cracks() {
    let c = [16.0, 16.0, 16.0];
    let field = blob(c, 9.0, 1.5);
    let blocks: Vec<_> = tiles(0, 32, 8, 0, None).into_iter().map(|b| (b, TransitionSides::NONE)).collect();
    let mesh = build(&field, &blocks, &ExtractOptions::default());
    assert_closed_outward(&mesh, c);

    let coarse: Vec<_> = tiles(0, 32, 4, 1, None).into_iter().map(|b| (b, TransitionSides::NONE)).collect();
    assert_closed_outward(&build(&field, &coarse, &ExtractOptions::default()), c);
}

/// A coarse block at `lod` covering `[0, 16)³` with transition cells on all six
/// faces, surrounded by finer blocks on every side.
fn surrounded_scene(lod: u8) -> Vec<(Block, TransitionSides)> {
    let coarse_cells = 16 >> lod;
    let fine_cells = 8 >> (lod - 1);
    let mut blocks = vec![(Block::new([0, 0, 0], coarse_cells, lod), TransitionSides::ALL)];
    for b in tiles(-8, 24, fine_cells, lod - 1, Some(([0; 3], [16; 3]))) {
        blocks.push((b, TransitionSides::NONE));
    }
    blocks
}

#[test]
fn transition_faces_edges_and_corners_are_crack_free() {
    let c = [8.0, 8.0, 8.0];
    // Radii that cross only faces, also the block's edges, and also its corners.
    for r in [9.3, 12.2, 14.6] {
        for lod in [1, 2] {
            for width in [0.25, 0.5, 0.8] {
                let options = ExtractOptions::default().with_transition_width(width);
                let mesh = build(&sphere(c, r), &surrounded_scene(lod), &options);
                if let Err(e) = check_closed(&mesh) {
                    panic!("r={r} lod={lod} width={width}: {e}");
                }
                assert!(signed_volume(&mesh) > 0.0);
                // Thin transition slabs can fold slightly where the projected
                // displacement slides vertices further than the slab is deep.
                if width >= 0.5 {
                    assert!(outward_fraction(&mesh, c) > 0.97, "r={r} lod={lod} width={width}");
                }
            }
        }
    }
}

#[test]
fn lumpy_surface_through_transitions_is_crack_free() {
    let c = [8.0, 8.0, 8.0];
    for (r, amp) in [(10.0, 2.0), (12.0, 2.5), (13.0, 1.0)] {
        for lod in [1, 2] {
            let mesh = build(&blob(c, r, amp), &surrounded_scene(lod), &ExtractOptions::default());
            if let Err(e) = check_closed(&mesh) {
                panic!("r={r} amp={amp} lod={lod}: {e}");
            }
            assert!(signed_volume(&mesh) > 0.0);
        }
    }
}

#[test]
fn single_transition_face() {
    // Coarse block with its -X face against four finer blocks.
    let c = [0.0, 8.0, 8.0];
    let mut blocks = vec![(Block::new([0, 0, 0], 8, 1), TransitionSides::from(Face::NegX))];
    for b in tiles(-8, 16, 8, 0, None) {
        if b.origin[0] == -8 && b.origin[1] >= 0 && b.origin[2] >= 0 {
            blocks.push((b, TransitionSides::NONE));
        }
    }
    assert_eq!(blocks.len(), 5);
    let mesh = build(&blob(c, 5.0, 1.0), &blocks, &ExtractOptions::default());
    assert_closed_outward(&mesh, c);
}

#[test]
fn surface_on_exact_voxel_corners() {
    // Integer-valued densities put many vertices exactly on lattice points,
    // exercising the corner-reuse paths of both cell types.
    let c = [8.0, 8.0, 8.0];
    let base = sphere(c, 11.0);
    let field = move |x, y, z| base(x, y, z).round();
    for lod in [1, 2] {
        let mesh = build(&field, &surrounded_scene(lod), &ExtractOptions::default());
        // The surface may run exactly along a block edge here, where sheets touch.
        if let Err(e) = check_watertight(&mesh) {
            panic!("lod={lod}: {e}");
        }
        assert!(signed_volume(&mesh) > 0.0);
    }

    // Within a single block every corner vertex is shared, never duplicated.
    let mesh =
        extract(&field, &Block::new([-8, -8, -8], 32, 0), TransitionSides::NONE, &ExtractOptions::default()).unwrap();
    let unique: HashSet<_> = mesh.positions.iter().map(|p| p.map(f32::to_bits)).collect();
    assert_eq!(unique.len(), mesh.vertex_count());
    check_closed(&mesh).unwrap();
}

#[test]
fn transition_mesh_vertices_are_not_duplicated() {
    let c = [8.0, 8.0, 8.0];
    let field = sphere(c, 12.2);
    let mesh = extract(&field, &Block::new([0, 0, 0], 8, 1), TransitionSides::ALL, &ExtractOptions::default()).unwrap();
    let unique: HashSet<_> = mesh.positions.iter().map(|p| p.map(f32::to_bits)).collect();
    assert_eq!(unique.len(), mesh.vertex_count());
}

#[test]
fn isovalue_and_voxel_grid() {
    let c = [8.0, 8.0, 8.0];
    let grid = VoxelGrid::from_fn([-2, -2, -2], [21, 21, 21], sphere(c, 0.0));
    let options = ExtractOptions::default().with_isovalue(5.5);
    let mesh = extract(&grid, &Block::new([0, 0, 0], 16, 0), TransitionSides::NONE, &options).unwrap();
    assert_closed_outward(&mesh, c);
    let v = signed_volume(&mesh);
    let exact = 4.0 / 3.0 * std::f64::consts::PI * 5.5f64.powi(3);
    assert!((v - exact).abs() / exact < 0.05);
}

#[test]
fn errors() {
    let f = |_: i32, _: i32, _: i32| 1.0;
    let o = ExtractOptions::default();
    assert_eq!(extract(&f, &Block::new([0; 3], 0, 0), TransitionSides::NONE, &o), Err(ExtractError::EmptyBlock));
    assert_eq!(extract(&f, &Block::new([0; 3], 4, 0), Face::NegX.into(), &o), Err(ExtractError::TransitionAtLodZero));
    assert_eq!(
        extract(&f, &Block::new([0; 3], 4, 1), Face::NegX.into(), &o.with_transition_width(1.0)),
        Err(ExtractError::InvalidTransitionWidth)
    );
    assert_eq!(
        extract(&f, &Block::new([i32::MAX - 4, 0, 0], 4, 0), TransitionSides::NONE, &o),
        Err(ExtractError::CoordinateOverflow)
    );
    assert!(extract(&f, &Block::new([0; 3], 4, 1), TransitionSides::ALL, &o).unwrap().is_empty());
}
