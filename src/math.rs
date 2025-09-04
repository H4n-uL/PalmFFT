use crate::{Complex, ComplexExt};
use core::f64::consts::PI;

pub fn my_sincosm1pi(a: f64) -> Complex {
    let mut s = a * a;
    let mut r = -1.0369917389758117e-4;
    r *= s; r += 1.9294935641298806e-3;
    r *= s; r += -2.5806887942825395e-2;
    r *= s; r += 2.3533063028328211e-1;
    r *= s; r += -1.3352627688538006e+0;
    r *= s; r += 4.0587121264167623e+0;
    r *= s; r += -4.9348022005446790e+0;
    let c = r * s;
    let mut r2 = 4.6151442520157035e-4;
    r2 *= s; r2 += -7.3700183130883555e-3;
    r2 *= s; r2 += 8.2145868949323936e-2;
    r2 *= s; r2 += -5.9926452893214921e-1;
    r2 *= s; r2 += 2.5501640398732688e+0;
    r2 *= s; r2 += -5.1677127800499516e+0;
    s = a * PI + (r2 * s * a);
    return Complex::new(c, s);
}

fn calc_first_octant(den: usize, res: &mut [Complex]) {
    let n = (den + 4) >> 3;
    if n == 0 { return; }
    res[0] = Complex::new(1.0, 0.0);
    if n == 1 { return; }

    let l1 = libm::sqrt(n as f64) as usize;
    for i in 1..l1 {
        res[i] = my_sincosm1pi((2.0 * (i as f64)) / (den as f64));
    }

    let mut start = l1;
    while start < n {
        let z = my_sincosm1pi((2.0 * (start as f64)) / (den as f64));
        res[start] = z + 1.0;

        let mut end = l1;
        if start + end > n { end = n - start; }
        for i in 1..end {
            res[start + i] = ((z * res[i] + z) + res[i]) + 1.0;
        }
        start += l1;
    }

    for i in 1..l1 { res[i] += 1.0; }
}

fn calc_first_quadrant(n: usize, res: &mut [Complex]) {
    let (head, p) = res.split_at_mut(n / 2);
    calc_first_octant(n << 1, p);
    let ndone = (n + 2) >> 2;
    let mut i = 0;
    let mut idx1 = 0;
    let mut idx2 = ndone - 1;

    while i + 1 < ndone {
        head[idx1] = p[i];
        head[idx2] = Complex::new(p[i + 1].im, p[i + 1].re);
        i += 2; idx1 += 1; idx2 -= 1;
    }
    if i != ndone {
        head[idx1] = p[i];
    }
}

fn calc_first_half(n: usize, res: &mut [Complex]) {
    let ndone = (n + 1) >> 1;
    let half = n >> 1;
    calc_first_octant(n << 2, res);
    res.rotate_right(half);

    let mut i4: isize = 0;
    let in_val = n as isize;
    let mut i = 0;

    while i4 <= in_val - i4 && i < ndone {
        res[i] = res[i4 as usize + half];
        i += 1; i4 += 4;
    }

    while i4 - in_val <= 0 && i < ndone {
        res[i] = res[(in_val - i4) as usize + half].rotm90().conj();
        i += 1; i4 += 4;
    }

    while i4 <= 3 * in_val - i4 && i < ndone {
        res[i] = res[(i4 - in_val) as usize + half].rot90();
        i += 1; i4 += 4;
    }

    while i < ndone {
        res[i] = -res[(2 * in_val - i4) as usize + half].conj();
        i += 1; i4 += 4;
    }
}

fn fill_first_quadrant(n: usize, res: &mut [Complex]) {
    let hsqt2 = 0.707106781186547524400844362104849_f64;
    let quart = n >> 2;
    let eighth = n >> 3;
    if (n & 7) == 0 {
        res[eighth] = Complex::new(hsqt2, hsqt2);
    }

    for (i, j) in (1..=eighth).zip((0..quart).rev()) {
        res[j] = Complex::new(res[i].im, res[i].re);
    }
}

fn fill_first_half(n: usize, res: &mut [Complex]) {
    let half = n >> 1;
    let quart = n >> 2;
    if (n & 3) == 0 {
        for i in 0..quart {
            res[quart + i] = Complex::new(-res[i].im, res[i].re);
        }
    }
    else {
        for (i, j) in (1..=quart).zip((0..half).rev()) {
            res[j] = Complex::new(-res[i].re, res[i].im);
        }
    }
}

fn fill_second_half(n: usize, res: &mut [Complex]) {
    let half = n >> 1;
    if (n & 1) == 0 {
        for i in 0..half {
            res[i + half] = -res[i];
        }
    }
    else {
        for (i, j) in (1..=half).zip((half..n).rev()) {
            res[j] = res[i].conj()
        }
    }
}

fn sincos_2pibyn_half(n: usize, res: &mut [Complex]) {
    if (n & 3) == 0 {
        calc_first_octant(n, res);
        fill_first_quadrant(n, res);
        fill_first_half(n, res);
    }
    else if (n & 1) == 0 {
        calc_first_quadrant(n, res);
        fill_first_half(n, res);
    }
    else { calc_first_half(n, res); }
}

pub fn sincos_2pibyn(n: usize, res: &mut [Complex]) {
    sincos_2pibyn_half(n, res);
    fill_second_half(n, res);
}

pub fn largest_prime_factor(mut n: usize) -> usize {
    let mut max_prime = 1;
    while n & 1 == 0 {
        max_prime = 2;
        n >>= 1;
    }
    let mut factor = 3;
    while factor * factor <= n {
        while n % factor == 0 {
            max_prime = factor;
            n /= factor;
        }
        factor += 2;
    }
    if n > 2 { max_prime = n; }
    return max_prime;
}

pub fn cost_guess(mut n: usize) -> f64 {
    const LFP: f64 = 1.1;
    let ni = n;
    let mut result = 0.0;
    while n & 1 == 0 {
        result += 2.0;
        n >>= 1;
    }
    let mut limit = libm::sqrt(n as f64 + 0.01) as usize;
    let mut x = 3; while x <= limit {
        while n % x == 0 {
            result += if x <= 5 { x as f64 } else { LFP * x as f64 };
            n /= x; limit = libm::sqrt(n as f64 + 0.01) as usize;
        }
    x += 2; }
    if n > 1 { result += if n <= 5 { n as f64 } else { LFP * n as f64 }; }
    return result * ni as f64;
}

pub fn good_size(n: usize) -> usize {
    if n <= 6 { return n; }
    let mut bestfac = 2 * n;

    let mut f2 = 1; while f2 < bestfac {
        let mut f23 = f2; while f23 < bestfac {
            let mut f235 = f23; while f235 < bestfac {
                let mut f2357 = f235; while f2357 < bestfac {
                    let mut f235711 = f2357; while f235711 < bestfac {
                        if f235711 >= n { bestfac = f235711; }
                    f235711 *= 11; }
                f2357 *= 7; }
            f235 *= 5; }
        f23 *= 3; }
    f2 *= 2; }

    return bestfac;
}