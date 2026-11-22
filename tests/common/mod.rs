#![allow(dead_code)]

use std::collections::HashMap;
use transvoxel::Mesh;

pub fn sphere(c: [f32; 3], r: f32) -> impl Fn(i32, i32, i32) -> f32 {
    move |x, y, z| {
        let d = [x as f32 - c[0], y as f32 - c[1], z as f32 - c[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - r
    }
}

/// A lumpy sphere: radius `r` with bumps of amplitude `amp`.
pub fn blob(c: [f32; 3], r: f32, amp: f32) -> impl Fn(i32, i32, i32) -> f32 {
    move |x, y, z| {
        let d = [x as f32 - c[0], y as f32 - c[1], z as f32 - c[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let bumps = (d[0] * 0.9).sin() * (d[1] * 0.7 + 1.0).cos() + (d[2] * 1.1 + 0.3).sin() * 0.5;
        len - r + amp * bumps
    }
}

fn key(p: [f32; 3]) -> [i64; 3] {
    [
        (p[0] as f64 * 65536.0).round() as i64,
        (p[1] as f64 * 65536.0).round() as i64,
        (p[2] as f64 * 65536.0).round() as i64,
    ]
}

/// Checks that the mesh is a closed, consistently oriented surface when vertices
/// are identified by position: every directed edge appears exactly once and its
/// reverse appears exactly once.
pub fn check_closed(mesh: &Mesh) -> Result<(), String> {
    check_edges(mesh, true)
}

/// Weaker than [`check_closed`]: every directed edge is matched by as many
/// reverse edges, so there are no holes, but sheets may touch along an edge.
pub fn check_watertight(mesh: &Mesh) -> Result<(), String> {
    check_edges(mesh, false)
}

fn check_edges(mesh: &Mesh, manifold: bool) -> Result<(), String> {
    let mut ids: HashMap<[i64; 3], usize> = HashMap::new();
    let canon: Vec<usize> = mesh
        .positions
        .iter()
        .map(|&p| {
            let n = ids.len();
            *ids.entry(key(p)).or_insert(n)
        })
        .collect();
    let mut edges: HashMap<(usize, usize), u32> = HashMap::new();
    let mut degenerate = 0;
    for t in mesh.triangles() {
        let c = [canon[t[0] as usize], canon[t[1] as usize], canon[t[2] as usize]];
        if c[0] == c[1] || c[1] == c[2] || c[0] == c[2] {
            degenerate += 1;
            continue;
        }
        for k in 0..3 {
            *edges.entry((c[k], c[(k + 1) % 3])).or_default() += 1;
        }
    }
    // Zero-area triangles (two distinct vertices at one position) occur where the
    // projected transition displacement vanishes. Dropping them removes a pair of
    // opposite edges, so they never hide a crack.
    let _ = degenerate;
    let mut bad = Vec::new();
    for (&(a, b), &count) in &edges {
        let back = edges.get(&(b, a)).copied().unwrap_or(0);
        if (manifold && (count != 1 || back != 1)) || count != back {
            bad.push((a, b, count, back));
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        let pos: Vec<[i64; 3]> = {
            let mut v = vec![[0; 3]; ids.len()];
            for (k, i) in &ids {
                v[*i] = *k;
            }
            v
        };
        let sample: Vec<String> = bad
            .iter()
            .take(6)
            .map(|&(a, b, c, r)| {
                let f = |k: [i64; 3]| [k[0] as f64 / 65536.0, k[1] as f64 / 65536.0, k[2] as f64 / 65536.0];
                format!("{:?}->{:?} fwd={c} back={r}", f(pos[a]), f(pos[b]))
            })
            .collect();
        Err(format!("{} bad directed edges of {}; e.g. {}", bad.len(), edges.len(), sample.join("; ")))
    }
}

/// Signed volume enclosed by the mesh (positive for outward-facing triangles).
pub fn signed_volume(mesh: &Mesh) -> f64 {
    let mut v = 0.0;
    for t in mesh.triangles() {
        let p: Vec<[f64; 3]> = t
            .iter()
            .map(|&i| {
                let q = mesh.positions[i as usize];
                [q[0] as f64, q[1] as f64, q[2] as f64]
            })
            .collect();
        let c = [
            p[1][1] * p[2][2] - p[1][2] * p[2][1],
            p[1][2] * p[2][0] - p[1][0] * p[2][2],
            p[1][0] * p[2][1] - p[1][1] * p[2][0],
        ];
        v += p[0][0] * c[0] + p[0][1] * c[1] + p[0][2] * c[2];
    }
    v / 6.0
}

/// Geometric normal of every triangle compared with the outward direction from
/// `center`: returns the fraction of triangles facing outward.
pub fn outward_fraction(mesh: &Mesh, center: [f32; 3]) -> f64 {
    let mut good = 0usize;
    let mut total = 0usize;
    for t in mesh.triangles() {
        let p = t.map(|i| mesh.positions[i as usize]);
        let e1 = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
        let e2 = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
        let n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        let c = [
            (p[0][0] + p[1][0] + p[2][0]) / 3.0 - center[0],
            (p[0][1] + p[1][1] + p[2][1]) / 3.0 - center[1],
            (p[0][2] + p[1][2] + p[2][2]) / 3.0 - center[2],
        ];
        let d = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];
        if n.iter().map(|v| v * v).sum::<f32>() < 1e-12 {
            continue;
        }
        total += 1;
        if d > 0.0 {
            good += 1;
        }
    }
    good as f64 / total.max(1) as f64
}
