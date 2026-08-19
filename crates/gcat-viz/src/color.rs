//! Turning counts into pixels.
//!
//! Writes RGBA straight into a caller-supplied buffer. In the browser that
//! buffer lives in WebAssembly memory and JavaScript wraps it in an `ImageData`
//! without copying, so a megapixel image costs one pass and no serialisation.

/// How counts map to colour.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub enum Palette {
    /// Dark ground, warm highlights. Reads well on a dark page.
    #[default]
    Ember,
    /// Light ground, cool highlights. Reads well on a light page.
    Ice,
    /// Plain grey ramp.
    Mono,
}

impl Palette {
    /// Map a normalised intensity in `0..=1` to RGB.
    fn sample(self, t: f32) -> (u8, u8, u8) {
        let t = t.clamp(0.0, 1.0);
        match self {
            // Piecewise-linear ramps. Three stops are enough for a histogram
            // and avoid the banding a two-stop lerp shows in the midtones.
            //
            // Ember starts at true black so empty cells vanish into a black
            // page rather than sitting on a faintly blue square, then climbs
            // through orange to gold.
            Palette::Ember => ramp(t, (0, 0, 0), (214, 62, 0), (255, 192, 20)),
            Palette::Ice => ramp(t, (255, 255, 255), (26, 168, 255), (0, 40, 90)),
            Palette::Mono => ramp(t, (16, 16, 16), (128, 128, 128), (245, 245, 245)),
        }
    }
}

fn ramp(t: f32, a: (u8, u8, u8), b: (u8, u8, u8), c: (u8, u8, u8)) -> (u8, u8, u8) {
    let (from, to, local) = if t < 0.5 {
        (a, b, t * 2.0)
    } else {
        (b, c, (t - 0.5) * 2.0)
    };
    let lerp = |x: u8, y: u8| (x as f32 + (y as f32 - x as f32) * local) as u8;
    (lerp(from.0, to.0), lerp(from.1, to.1), lerp(from.2, to.2))
}

/// Write `counts` into `rgba` as a colourised image.
///
/// Intensity is scaled by `count^gamma / peak^gamma`. Genomic count
/// distributions are heavy-tailed — a handful of repeat k-mers can be orders of
/// magnitude above the median — so a linear ramp renders as a black field with
/// three bright dots. A gamma below 1 pulls the bulk of the distribution into
/// the visible range.
///
/// `rgba` must hold `4 * counts.len()` bytes; extra capacity is ignored.
pub fn colorize(counts: &[u32], peak: u32, palette: Palette, gamma: f32, rgba: &mut [u8]) {
    let usable = counts.len().min(rgba.len() / 4);
    if peak == 0 {
        let (r, g, b) = palette.sample(0.0);
        for i in 0..usable {
            rgba[i * 4] = r;
            rgba[i * 4 + 1] = g;
            rgba[i * 4 + 2] = b;
            rgba[i * 4 + 3] = 255;
        }
        return;
    }

    let inv_peak = 1.0 / peak as f32;
    for i in 0..usable {
        let c = counts[i];
        let t = if c == 0 {
            0.0
        } else {
            powf(c as f32 * inv_peak, gamma)
        };
        let (r, g, b) = palette.sample(t);
        rgba[i * 4] = r;
        rgba[i * 4 + 1] = g;
        rgba[i * 4 + 2] = b;
        rgba[i * 4 + 3] = 255;
    }
}

/// `x^p` for `x` in `0..=1`, without `std`.
fn powf(x: f32, p: f32) -> f32 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    exp(p * ln(x))
}

fn ln(x: f32) -> f32 {
    let bits = x.to_bits();
    let exp_part = ((bits >> 23) & 0xff) as i32 - 127;
    let mantissa = f32::from_bits((bits & 0x007f_ffff) | 0x3f80_0000);
    let z = (mantissa - 1.0) / (mantissa + 1.0);
    let z2 = z * z;
    let mut series = 0.0f32;
    let mut term = z;
    let mut k = 1.0f32;
    for _ in 0..12 {
        series += term / k;
        term *= z2;
        k += 2.0;
    }
    core::f32::consts::LN_2 * exp_part as f32 + 2.0 * series
}

fn exp(x: f32) -> f32 {
    if x < -30.0 {
        return 0.0;
    }
    if x > 30.0 {
        return f32::INFINITY;
    }
    // Range-reduce to |r| <= ln(2)/2, where the Taylor series converges
    // quickly, then scale back by a power of two.
    let n = (x * core::f32::consts::LOG2_E + 0.5 * if x < 0.0 { -1.0 } else { 1.0 }) as i32;
    let r = x - n as f32 * core::f32::consts::LN_2;
    let mut term = 1.0f32;
    let mut sum = 1.0f32;
    for i in 1..12 {
        term *= r / i as f32;
        sum += term;
    }
    sum * pow2i(n)
}

fn pow2i(n: i32) -> f32 {
    let clamped = n.clamp(-126, 127);
    f32::from_bits(((clamped + 127) as u32) << 23)
}
