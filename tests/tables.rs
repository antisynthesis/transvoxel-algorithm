//! Structural invariants of the lookup tables.

use transvoxel::tables::*;

#[test]
fn regular_tables_are_consistent() {
    for case in 0..256usize {
        let data = &REGULAR_CELL_DATA[REGULAR_CELL_CLASS[case] as usize];
        let verts = &REGULAR_VERTEX_DATA[case];
        let count = data.vertex_count();
        assert!(verts[..count].iter().all(|&v| v != 0), "case {case}");
        assert!(verts[count..].iter().all(|&v| v == 0), "case {case}");
        assert!(data.indices().iter().all(|&i| (i as usize) < count), "case {case}");
        let mut used = vec![false; count];
        data.indices().iter().for_each(|&i| used[i as usize] = true);
        assert!(used.iter().all(|&u| u), "case {case}: unused vertex");

        for &code in &verts[..count] {
            let v0 = ((code >> 4) & 0xF) as usize;
            let v1 = (code & 0xF) as usize;
            assert!(v0 < v1 && v1 < 8);
            assert_eq!((v0 ^ v1).count_ones(), 1, "case {case}: not an edge");
            let inside = |c: usize| case >> c & 1 == 1;
            assert_ne!(inside(v0), inside(v1), "case {case}: edge does not cross");

            // The reuse code must name the owning cell and the slot of that edge in it.
            let dir = (code >> 12) as usize;
            let slot = ((code >> 8) & 0xF) as usize;
            let owned = match slot {
                1 => (5, 7),
                2 => (6, 7),
                3 => (3, 7),
                _ => panic!("case {case}: bad slot {slot}"),
            };
            if dir == 8 {
                assert_eq!((v0, v1), owned, "case {case}");
            } else {
                assert!((1..8).contains(&dir));
                assert_eq!(v0 & dir, 0);
                assert_eq!((v0 | dir, v1 | dir), owned, "case {case}");
            }
        }
    }
}

#[test]
fn every_class_is_used() {
    let mut regular = [false; 16];
    REGULAR_CELL_CLASS.iter().for_each(|&c| regular[c as usize] = true);
    assert!(regular.iter().all(|&u| u));
    let mut transition = [false; 56];
    TRANSITION_CELL_CLASS.iter().for_each(|&c| transition[(c & 0x7F) as usize] = true);
    assert!(transition.iter().all(|&u| u));
}

#[test]
fn transition_tables_are_consistent() {
    for case in 0..512usize {
        let class = TRANSITION_CELL_CLASS[case];
        let data = &TRANSITION_CELL_DATA[(class & 0x7F) as usize];
        let verts = &TRANSITION_VERTEX_DATA[case];
        let count = data.vertex_count();
        assert!(count <= 12);
        assert!(verts[..count].iter().all(|&v| v != 0), "case {case}");
        assert!(verts[count..].iter().all(|&v| v == 0), "case {case}");
        assert!(data.indices().iter().all(|&i| (i as usize) < count), "case {case}");

        let mut inside = [false; 13];
        for (bit, &s) in TRANSITION_CASE_SAMPLES.iter().enumerate() {
            inside[s] = case >> bit & 1 == 1;
        }
        inside[9] = inside[0];
        inside[10] = inside[2];
        inside[11] = inside[6];
        inside[12] = inside[8];

        for &code in &verts[..count] {
            let s0 = ((code >> 4) & 0xF) as usize;
            let s1 = (code & 0xF) as usize;
            assert!(s0 < 13 && s1 < 13);
            assert_eq!(s0 >= 9, s1 >= 9, "case {case}: edge spans both faces");
            assert_ne!(inside[s0], inside[s1], "case {case}: edge does not cross");
            let dir = code >> 12;
            assert!(matches!(dir, 1 | 2 | 4 | 8), "case {case}: dir {dir}");
        }
    }
    assert_eq!(TRANSITION_CELL_CLASS[0], 0);
    assert_eq!(TRANSITION_CELL_CLASS[511] & 0x7F, 0);
}

#[test]
fn transition_corner_data_matches_ownership() {
    // Slot 0 = sample 8, 1 = sample 7, 2 = sample 5, 7 = sample C.
    assert_eq!(TRANSITION_CORNER_DATA[8], 0x80);
    assert_eq!(TRANSITION_CORNER_DATA[7], 0x81);
    assert_eq!(TRANSITION_CORNER_DATA[5], 0x82);
    assert_eq!(TRANSITION_CORNER_DATA[12], 0x87);
    assert_eq!(TRANSITION_CORNER_DATA[4] >> 4, 4);
}
