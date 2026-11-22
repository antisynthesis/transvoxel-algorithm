//! Builds a two-level terrain and writes it as a Wavefront OBJ file.
//!
//! A 3×3×3 grid of LOD-1 blocks covers voxels `0..96`; the centre block is
//! replaced by eight LOD-0 blocks. Each coarse block that touches the centre
//! gets a transition face there, so the result has no cracks.
//!
//! ```text
//! cargo run --release --example terrain -- terrain.obj
//! ```

use std::fmt::Write as _;
use transvoxel::{extract_into, Block, ExtractOptions, Face, Mesh, TransitionSides};

fn terrain(x: i32, y: i32, z: i32) -> f32 {
    let (x, y, z) = (x as f32, y as f32, z as f32);
    let height = 48.0 + 9.0 * (x * 0.07).sin() * (z * 0.05).cos() + 4.0 * (x * 0.13 + z * 0.11).sin();
    let cave = 12.0 - ((x - 48.0).powi(2) + (y - 44.0).powi(2) * 2.0 + (z - 48.0).powi(2)).sqrt();
    // Solid (negative) below the height field, with a cave carved out.
    (y - height).max(cave)
}

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| "terrain.obj".into());
    let options = ExtractOptions::default();
    let mut mesh = Mesh::new();

    let centre = [32, 32, 32];
    for bz in 0..3 {
        for by in 0..3 {
            for bx in 0..3 {
                let origin = [bx * 32, by * 32, bz * 32];
                if origin == centre {
                    for k in 0..8 {
                        let o = [32 + (k & 1) * 16, 32 + ((k >> 1) & 1) * 16, 32 + ((k >> 2) & 1) * 16];
                        extract_into(&terrain, &Block::new(o, 16, 0), TransitionSides::NONE, &options, &mut mesh)
                            .unwrap();
                    }
                    continue;
                }
                // Faces of this coarse block that touch the fine centre block.
                let mut sides = TransitionSides::NONE;
                let b = [bx, by, bz];
                for face in Face::ALL {
                    let mut n = b;
                    n[face.axis()] += if face.is_positive() { 1 } else { -1 };
                    if n == [1, 1, 1] {
                        sides |= face;
                    }
                }
                extract_into(&terrain, &Block::new(origin, 16, 1), sides, &options, &mut mesh).unwrap();
            }
        }
    }

    let mut obj = String::new();
    for p in &mesh.positions {
        writeln!(obj, "v {} {} {}", p[0], p[1], p[2]).unwrap();
    }
    for n in &mesh.normals {
        writeln!(obj, "vn {} {} {}", n[0], n[1], n[2]).unwrap();
    }
    for t in mesh.triangles() {
        let [a, b, c] = t.map(|i| i + 1);
        writeln!(obj, "f {a}//{a} {b}//{b} {c}//{c}").unwrap();
    }
    std::fs::write(&path, obj).expect("failed to write OBJ");
    println!("wrote {} vertices, {} triangles to {path}", mesh.vertex_count(), mesh.triangle_count());
}
