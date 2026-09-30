//! Gaussian quadrature rules for finite element integration.
//!
//! Implements 2- and 3-point 1D Gauss-Legendre rules and their 2x2 and 3x3
//! tensor-product forms. Citation gap: no literature source is cited for these
//! tables (the MITC3 Hammer rule lives in `elements/mitc3.rs`).
//!
//! All rules return statically sized point/weight arrays to avoid heap
//! allocation in the hot integration path.

// ============================================================================
// Element-specific rules (exported as pub const arrays for use in element files)
// ============================================================================

// --- 1D Gauss rules on [-1,1] ---

/// 1D Gauss points for 2-point rule: ±1/√3.
pub const GAUSS2_PTS: [f64; 2] = [-GP, GP];
/// 1D Gauss weights for 2-point rule: [1.0, 1.0].
pub const GAUSS2_W: [f64; 2] = [1.0, 1.0];

/// 1D Gauss points for 3-point rule.
pub const GAUSS3_PTS: [f64; 3] = [-SQRT35, 0.0, SQRT35];
/// 1D Gauss weights for 3-point rule.
pub const GAUSS3_W: [f64; 3] = [W5_9, W8_9, W5_9];

// ============================================================================
// Static 2D Gauss point/weight arrays (for ReferenceElement2D trait impls)
// ============================================================================

/// 2×2 Gauss rule on [−1,1]²: 4 points.
pub static GAUSS2X2_PTS: [[f64; 2]; 4] = [
    [-0.577_350_269_189_625_8, -0.577_350_269_189_625_8],
    [ 0.577_350_269_189_625_8, -0.577_350_269_189_625_8],
    [ 0.577_350_269_189_625_8,  0.577_350_269_189_625_8],
    [-0.577_350_269_189_625_8,  0.577_350_269_189_625_8],
];
/// 2×2 Gauss weights (all 1.0).
pub static GAUSS2X2_W: [f64; 4] = [1.0, 1.0, 1.0, 1.0];

/// 3×3 Gauss rule on [−1,1]²: 9 points.
pub static GAUSS3X3_PTS: [[f64; 2]; 9] = [
    [-0.774_596_669_241_483_4, -0.774_596_669_241_483_4],
    [ 0.0,                     -0.774_596_669_241_483_4],
    [ 0.774_596_669_241_483_4, -0.774_596_669_241_483_4],
    [-0.774_596_669_241_483_4,  0.0                    ],
    [ 0.0,                      0.0                    ],
    [ 0.774_596_669_241_483_4,  0.0                    ],
    [-0.774_596_669_241_483_4,  0.774_596_669_241_483_4],
    [ 0.0,                      0.774_596_669_241_483_4],
    [ 0.774_596_669_241_483_4,  0.774_596_669_241_483_4],
];
/// 3×3 Gauss weights (W5/9 × W5/9, etc.).
pub static GAUSS3X3_W: [f64; 9] = [
    0.555_555_555_555_555_6 * 0.555_555_555_555_555_6,
    0.888_888_888_888_888_9 * 0.555_555_555_555_555_6,
    0.555_555_555_555_555_6 * 0.555_555_555_555_555_6,
    0.555_555_555_555_555_6 * 0.888_888_888_888_888_9,
    0.888_888_888_888_888_9 * 0.888_888_888_888_888_9,
    0.555_555_555_555_555_6 * 0.888_888_888_888_888_9,
    0.555_555_555_555_555_6 * 0.555_555_555_555_555_6,
    0.888_888_888_888_888_9 * 0.555_555_555_555_555_6,
    0.555_555_555_555_555_6 * 0.555_555_555_555_555_6,
];

// ============================================================================
// Internal constants
// ============================================================================

/// 1/√3 (2-pt Gauss point on [-1,1]).
const GP: f64 = 0.577_350_269_189_625_8;
/// √(3/5) (3-pt Gauss point on [-1,1]).
const SQRT35: f64 = 0.774_596_669_241_483_4;
/// 5/9 (3-pt Gauss weight for end points).
const W5_9: f64 = 0.555_555_555_555_555_6;
/// 8/9 (3-pt Gauss weight for midpoint).
const W8_9: f64 = 0.888_888_888_888_888_9;
