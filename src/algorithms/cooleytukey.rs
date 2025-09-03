use crate::{
    math::sincos_2pibyn, Complex, ComplexExt
};

fn pmc(a: Complex, b: Complex) -> (Complex, Complex) { (a + b, a - b) }

fn pass2(ido: usize, l1: usize, cc: &[Complex], ch: &mut [Complex], wa: &[Complex], sign: i32) {
    let cdim = 2;
    if ido == 1 {
        for k in 0..l1 {
            (ch[k], ch[k + l1]) = pmc(cc[k * cdim], cc[k * cdim + 1]);
        }
    }
    else {
        for k in 0..l1 {
            (ch[k * ido], ch[(k + l1) * ido]) = pmc(cc[k * cdim * ido], cc[k * cdim * ido + 1 * ido]);

            for i in 1..ido {
                ch[i + k * ido] = cc[i + k * cdim * ido] + cc[i + (k * cdim + 1) * ido];
                let twiddle = if sign > 0 { wa[i - 1] } else { wa[i - 1].conj() };
                ch[i + (k + l1) * ido] = twiddle * (cc[i + k * cdim * ido] - cc[i + (k * cdim + 1) * ido]);
            }
        }
    }
}

fn pass3(ido: usize, l1: usize, cc: &[Complex], ch: &mut [Complex], wa: &[Complex], sign: i32) {
    let cdim = 3;
    const TW1R: f64 = -0.5;
    let tw1i = (sign as f64) * 0.86602540378443864676;

    if ido == 1 {
        for k in 0..l1 {
            let t0 = cc[k * cdim];
            let (t1, t2) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (2 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t0 + t1;

            let ca = Complex::new(t0.re + TW1R * t1.re, t0.im + TW1R * t1.im);
            let cb = Complex::new(-tw1i * t2.im, tw1i * t2.re);
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 2)]) = pmc(ca, cb);
        }
    }
    else {
        for k in 0..l1 {
            let t0 = cc[0 + ido * (0 + k * cdim)];
            let (t1, t2) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (2 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t0 + t1;

            let ca = Complex::new(t0.re + TW1R * t1.re, t0.im + TW1R * t1.im);
            let cb = Complex::new(-tw1i * t2.im, tw1i * t2.re);
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 2)]) = pmc(ca, cb);

            for i in 1..ido {
                let t1_in = cc[i + ido * (1 + k * cdim)];
                let t2_in = cc[i + ido * (2 + k * cdim)];
                let t0 = cc[i + ido * (0 + k * cdim)];
                let (t1, t2) = pmc(t1_in, t2_in);

                ch[i + ido * (k + l1 * 0)] = t0 + t1;

                let ca = Complex::new(t0.re + TW1R * t1.re, t0.im + TW1R * t1.im);
                let cb = Complex::new(-tw1i * t2.im, tw1i * t2.re);
                let (da, db) = pmc(ca, cb);

                if sign < 0 {
                    ch[i + ido * (k + l1 * 1)] = wa[i - 1 + 0 * (ido - 1)].conj() * da;
                    ch[i + ido * (k + l1 * 2)] = wa[i - 1 + 1 * (ido - 1)].conj() * db;
                }
                else {
                    ch[i + ido * (k + l1 * 1)] = wa[i - 1 + 0 * (ido - 1)] * da;
                    ch[i + ido * (k + l1 * 2)] = wa[i - 1 + 1 * (ido - 1)] * db;
                }
            }
        }
    }
}

fn pass4(ido: usize, l1: usize, cc: &[Complex], ch: &mut [Complex], wa: &[Complex], sign: i32) {
    let cdim = 4;

    if ido == 1 {
        for k in 0..l1 {
            let (t2, t1) = pmc(cc[k * cdim], cc[k * cdim + 2]);
            let (t3, mut t4) = pmc(cc[k * cdim + 1], cc[k * cdim + 3]);

            t4 = if sign < 0 { t4.rotm90() } else { t4.rot90() };

            (ch[0 + ido * (k + l1 * 0)], ch[0 + ido * (k + l1 * 2)]) = pmc(t2, t3);
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 3)]) = pmc(t1, t4);
        }
    }
    else {
        for k in 0..l1 {
            let (t2, t1) = pmc(cc[k * cdim * ido], cc[k * cdim * ido + 2 * ido]);
            let (t3, mut t4) = pmc(cc[k * cdim * ido + ido], cc[k * cdim * ido + 3 * ido]);

            t4 = if sign < 0 { t4.rotm90() } else { t4.rot90() };

            (ch[0 + ido * (k + l1 * 0)], ch[0 + ido * (k + l1 * 2)]) = pmc(t2, t3);
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 3)]) = pmc(t1, t4);

            for i in 1..ido {
                let (t2, t1) = pmc(cc[i + k * cdim * ido], cc[i + k * cdim * ido + 2 * ido]);
                let (t3, mut t4) = pmc(cc[i + k * cdim * ido + ido], cc[i + k * cdim * ido + 3 * ido]);

                t4 = if sign < 0 { t4.rotm90() } else { t4.rot90() };

                let (c0, c3) = pmc(t2, t3);
                let (c2, c4) = pmc(t1, t4);

                ch[i + ido * (k + l1 * 0)] = c0;
                if sign < 0 {
                    ch[i + ido * (k + l1 * 1)] = wa[i - 1 + 0 * (ido - 1)].conj() * c2;
                    ch[i + ido * (k + l1 * 2)] = wa[i - 1 + 1 * (ido - 1)].conj() * c3;
                    ch[i + ido * (k + l1 * 3)] = wa[i - 1 + 2 * (ido - 1)].conj() * c4;
                }
                else {
                    ch[i + ido * (k + l1 * 1)] = wa[i - 1 + 0 * (ido - 1)] * c2;
                    ch[i + ido * (k + l1 * 2)] = wa[i - 1 + 1 * (ido - 1)] * c3;
                    ch[i + ido * (k + l1 * 3)] = wa[i - 1 + 2 * (ido - 1)] * c4;
                }
            }
        }
    }
}

fn pass5(ido: usize, l1: usize, cc: &[Complex], ch: &mut [Complex], wa: &[Complex], sign: i32) {
    let cdim = 5;
    const TW1R: f64 = 0.3090169943749474241;
    const TW2R: f64 = -0.8090169943749474241;
    let tw1i = (sign as f64) * 0.95105651629515357212;
    let tw2i = (sign as f64) * 0.58778525229247312917;

    if ido == 1 {
        for k in 0..l1 {
            let t0 = cc[0 + ido * (0 + k * cdim)];
            let (t1, t4) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (4 + k * cdim)]);
            let (t2, t3) = pmc(cc[0 + ido * (2 + k * cdim)], cc[0 + ido * (3 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t0 + t1 + t2;

            let ca = Complex::new(
                t0.re + TW1R * t1.re + TW2R * t2.re,
                t0.im + TW1R * t1.im + TW2R * t2.im
            );
            let cb = Complex::new(
                -(tw1i * t4.im + tw2i * t3.im),
                tw1i * t4.re + tw2i * t3.re
            );
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 4)]) = pmc(ca, cb);

            let ca = Complex::new(
                t0.re + TW2R * t1.re + TW1R * t2.re,
                t0.im + TW2R * t1.im + TW1R * t2.im
            );
            let cb = Complex::new(
                -(tw2i * t4.im - tw1i * t3.im),
                tw2i * t4.re - tw1i * t3.re
            );
            (ch[0 + ido * (k + l1 * 2)], ch[0 + ido * (k + l1 * 3)]) = pmc(ca, cb);
        }
    }
    else {
        for k in 0..l1 {
            let t0 = cc[0 + ido * (0 + k * cdim)];
            let (t1, t4) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (4 + k * cdim)]);
            let (t2, t3) = pmc(cc[0 + ido * (2 + k * cdim)], cc[0 + ido * (3 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t0 + t1 + t2;

            let ca = Complex::new(
                t0.re + TW1R * t1.re + TW2R * t2.re,
                t0.im + TW1R * t1.im + TW2R * t2.im
            );
            let cb = Complex::new(
                -(tw1i * t4.im + tw2i * t3.im),
                tw1i * t4.re + tw2i * t3.re
            );
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 4)]) = pmc(ca, cb);

            let ca = Complex::new(
                t0.re + TW2R * t1.re + TW1R * t2.re,
                t0.im + TW2R * t1.im + TW1R * t2.im
            );
            let cb = Complex::new(
                -(tw2i * t4.im - tw1i * t3.im),
                tw2i * t4.re - tw1i * t3.re
            );
            (ch[0 + ido * (k + l1 * 2)], ch[0 + ido * (k + l1 * 3)]) = pmc(ca, cb);

            for i in 1..ido {
                let t0 = cc[i + ido * (0 + k * cdim)];
                let (t1, t4) = pmc(cc[i + ido * (1 + k * cdim)], cc[i + ido * (4 + k * cdim)]);
                let (t2, t3) = pmc(cc[i + ido * (2 + k * cdim)], cc[i + ido * (3 + k * cdim)]);

                ch[i + ido * (k + l1 * 0)] = t0 + t1 + t2;

                let ca = Complex::new(
                    t0.re + TW1R * t1.re + TW2R * t2.re,
                    t0.im + TW1R * t1.im + TW2R * t2.im
                );
                let cb = Complex::new(
                    -(tw1i * t4.im + tw2i * t3.im),
                    tw1i * t4.re + tw2i * t3.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 1)] = if sign < 0 { wa[i - 1 + 0 * (ido - 1)].conj() * da } else { wa[i - 1 + 0 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 4)] = if sign < 0 { wa[i - 1 + 3 * (ido - 1)].conj() * db } else { wa[i - 1 + 3 * (ido - 1)] * db };

                let ca = Complex::new(
                    t0.re + TW2R * t1.re + TW1R * t2.re,
                    t0.im + TW2R * t1.im + TW1R * t2.im
                );
                let cb = Complex::new(
                    -(tw2i * t4.im - tw1i * t3.im),
                    tw2i * t4.re - tw1i * t3.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 2)] = if sign < 0 { wa[i - 1 + 1 * (ido - 1)].conj() * da } else { wa[i - 1 + 1 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 3)] = if sign < 0 { wa[i - 1 + 2 * (ido - 1)].conj() * db } else { wa[i - 1 + 2 * (ido - 1)] * db };
            }
        }
    }
}

fn pass7(ido: usize, l1: usize, cc: &[Complex], ch: &mut [Complex], wa: &[Complex], sign: i32) {
    let cdim = 7;
    const TW1R: f64 = 0.623489801858733530525;
    const TW2R: f64 = -0.222520933956314404289;
    const TW3R: f64 = -0.9009688679024191262361;
    let tw1i = (sign as f64) * 0.7818314824680298087084;
    let tw2i = (sign as f64) * 0.9749279121818236070181;
    let tw3i = (sign as f64) * 0.4338837391175581204758;

    if ido == 1 {
        for k in 0..l1 {
            let t1 = cc[0 + ido * (0 + k * cdim)];
            let (t2, t7) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (6 + k * cdim)]);
            let (t3, t6) = pmc(cc[0 + ido * (2 + k * cdim)], cc[0 + ido * (5 + k * cdim)]);
            let (t4, t5) = pmc(cc[0 + ido * (3 + k * cdim)], cc[0 + ido * (4 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t1 + t2 + t3 + t4;

            let ca = Complex::new(
                t1.re + TW1R * t2.re + TW2R * t3.re + TW3R * t4.re,
                t1.im + TW1R * t2.im + TW2R * t3.im + TW3R * t4.im
            );
            let cb = Complex::new(
                -(tw1i * t7.im + tw2i * t6.im + tw3i * t5.im),
                tw1i * t7.re + tw2i * t6.re + tw3i * t5.re
            );
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 6)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW2R * t2.re + TW3R * t3.re + TW1R * t4.re,
                t1.im + TW2R * t2.im + TW3R * t3.im + TW1R * t4.im
            );
            let cb = Complex::new(
                -(tw2i * t7.im - tw3i * t6.im - tw1i * t5.im),
                tw2i * t7.re - tw3i * t6.re - tw1i * t5.re
            );
            (ch[0 + ido * (k + l1 * 2)], ch[0 + ido * (k + l1 * 5)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW3R * t2.re + TW1R * t3.re + TW2R * t4.re,
                t1.im + TW3R * t2.im + TW1R * t3.im + TW2R * t4.im
            );
            let cb = Complex::new(
                -(tw3i * t7.im - tw1i * t6.im + tw2i * t5.im),
                tw3i * t7.re - tw1i * t6.re + tw2i * t5.re
            );
            (ch[0 + ido * (k + l1 * 3)], ch[0 + ido * (k + l1 * 4)]) = pmc(ca, cb);
        }
    }
    else {
        for k in 0..l1 {
            let t1 = cc[0 + ido * (0 + k * cdim)];
            let (t2, t7) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (6 + k * cdim)]);
            let (t3, t6) = pmc(cc[0 + ido * (2 + k * cdim)], cc[0 + ido * (5 + k * cdim)]);
            let (t4, t5) = pmc(cc[0 + ido * (3 + k * cdim)], cc[0 + ido * (4 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t1 + t2 + t3 + t4;

            let ca = Complex::new(
                t1.re + TW1R * t2.re + TW2R * t3.re + TW3R * t4.re,
                t1.im + TW1R * t2.im + TW2R * t3.im + TW3R * t4.im
            );
            let cb = Complex::new(
                -(tw1i * t7.im + tw2i * t6.im + tw3i * t5.im),
                tw1i * t7.re + tw2i * t6.re + tw3i * t5.re
            );
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 6)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW2R * t2.re + TW3R * t3.re + TW1R * t4.re,
                t1.im + TW2R * t2.im + TW3R * t3.im + TW1R * t4.im
            );
            let cb = Complex::new(
                -(tw2i * t7.im - tw3i * t6.im - tw1i * t5.im),
                tw2i * t7.re - tw3i * t6.re - tw1i * t5.re
            );
            (ch[0 + ido * (k + l1 * 2)], ch[0 + ido * (k + l1 * 5)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW3R * t2.re + TW1R * t3.re + TW2R * t4.re,
                t1.im + TW3R * t2.im + TW1R * t3.im + TW2R * t4.im
            );
            let cb = Complex::new(
                -(tw3i * t7.im - tw1i * t6.im + tw2i * t5.im),
                tw3i * t7.re - tw1i * t6.re + tw2i * t5.re
            );
            (ch[0 + ido * (k + l1 * 3)], ch[0 + ido * (k + l1 * 4)]) = pmc(ca, cb);

            for i in 1..ido {
                let t1 = cc[i + ido * (0 + k * cdim)];
                let (t2, t7) = pmc(cc[i + ido * (1 + k * cdim)], cc[i + ido * (6 + k * cdim)]);
                let (t3, t6) = pmc(cc[i + ido * (2 + k * cdim)], cc[i + ido * (5 + k * cdim)]);
                let (t4, t5) = pmc(cc[i + ido * (3 + k * cdim)], cc[i + ido * (4 + k * cdim)]);

                ch[i + ido * (k + l1 * 0)] = t1 + t2 + t3 + t4;

                let ca = Complex::new(
                    t1.re + TW1R * t2.re + TW2R * t3.re + TW3R * t4.re,
                    t1.im + TW1R * t2.im + TW2R * t3.im + TW3R * t4.im
                );
                let cb = Complex::new(
                    -(tw1i * t7.im + tw2i * t6.im + tw3i * t5.im),
                    tw1i * t7.re + tw2i * t6.re + tw3i * t5.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 1)] = if sign < 0 { wa[i - 1 + 0 * (ido - 1)].conj() * da } else { wa[i - 1 + 0 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 6)] = if sign < 0 { wa[i - 1 + 5 * (ido - 1)].conj() * db } else { wa[i - 1 + 5 * (ido - 1)] * db };

                let ca = Complex::new(
                    t1.re + TW2R * t2.re + TW3R * t3.re + TW1R * t4.re,
                    t1.im + TW2R * t2.im + TW3R * t3.im + TW1R * t4.im
                );
                let cb = Complex::new(
                    -(tw2i * t7.im - tw3i * t6.im - tw1i * t5.im),
                    tw2i * t7.re - tw3i * t6.re - tw1i * t5.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 2)] = if sign < 0 { wa[i - 1 + 1 * (ido - 1)].conj() * da } else { wa[i - 1 + 1 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 5)] = if sign < 0 { wa[i - 1 + 4 * (ido - 1)].conj() * db } else { wa[i - 1 + 4 * (ido - 1)] * db };

                let ca = Complex::new(
                    t1.re + TW3R * t2.re + TW1R * t3.re + TW2R * t4.re,
                    t1.im + TW3R * t2.im + TW1R * t3.im + TW2R * t4.im
                );
                let cb = Complex::new(
                    -(tw3i * t7.im - tw1i * t6.im + tw2i * t5.im),
                    tw3i * t7.re - tw1i * t6.re + tw2i * t5.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 3)] = if sign < 0 { wa[i - 1 + 2 * (ido - 1)].conj() * da } else { wa[i - 1 + 2 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 4)] = if sign < 0 { wa[i - 1 + 3 * (ido - 1)].conj() * db } else { wa[i - 1 + 3 * (ido - 1)] * db };
            }
        }
    }
}

fn pass11(ido: usize, l1: usize, cc: &[Complex], ch: &mut [Complex], wa: &[Complex], sign: i32) {
    let cdim = 11;
    const TW1R: f64 = 0.8412535328311811688618;
    const TW2R: f64 = 0.4154150130018864255293;
    const TW3R: f64 = -0.1423148382732851404438;
    const TW4R: f64 = -0.6548607339452850640569;
    const TW5R: f64 = -0.9594929736144973898904;
    let tw1i = (sign as f64) * 0.5406408174555975821076;
    let tw2i = (sign as f64) * 0.9096319953545183714117;
    let tw3i = (sign as f64) * 0.9898214418809327323761;
    let tw4i = (sign as f64) * 0.755749574354258283774;
    let tw5i = (sign as f64) * 0.2817325568414296977114;

    if ido == 1 {
        for k in 0..l1 {
            let t1 = cc[0 + ido * (0 + k * cdim)];
            let (t2, t11) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (10 + k * cdim)]);
            let (t3, t10) = pmc(cc[0 + ido * (2 + k * cdim)], cc[0 + ido * (9 + k * cdim)]);
            let (t4, t9) = pmc(cc[0 + ido * (3 + k * cdim)], cc[0 + ido * (8 + k * cdim)]);
            let (t5, t8) = pmc(cc[0 + ido * (4 + k * cdim)], cc[0 + ido * (7 + k * cdim)]);
            let (t6, t7) = pmc(cc[0 + ido * (5 + k * cdim)], cc[0 + ido * (6 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t1 + t2 + t3 + t4 + t5 + t6;

            let ca = Complex::new(
                t1.re + TW1R * t2.re + TW2R * t3.re + TW3R * t4.re + TW4R * t5.re + TW5R * t6.re,
                t1.im + TW1R * t2.im + TW2R * t3.im + TW3R * t4.im + TW4R * t5.im + TW5R * t6.im
            );
            let cb = Complex::new(
                -(tw1i * t11.im + tw2i * t10.im + tw3i * t9.im + tw4i * t8.im + tw5i * t7.im),
                tw1i * t11.re + tw2i * t10.re + tw3i * t9.re + tw4i * t8.re + tw5i * t7.re
            );
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 10)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW2R * t2.re + TW4R * t3.re + TW5R * t4.re + TW3R * t5.re + TW1R * t6.re,
                t1.im + TW2R * t2.im + TW4R * t3.im + TW5R * t4.im + TW3R * t5.im + TW1R * t6.im
            );
            let cb = Complex::new(
                -(tw2i * t11.im + tw4i * t10.im - tw5i * t9.im - tw3i * t8.im - tw1i * t7.im),
                tw2i * t11.re + tw4i * t10.re - tw5i * t9.re - tw3i * t8.re - tw1i * t7.re
            );
            (ch[0 + ido * (k + l1 * 2)], ch[0 + ido * (k + l1 * 9)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW3R * t2.re + TW5R * t3.re + TW2R * t4.re + TW1R * t5.re + TW4R * t6.re,
                t1.im + TW3R * t2.im + TW5R * t3.im + TW2R * t4.im + TW1R * t5.im + TW4R * t6.im
            );
            let cb = Complex::new(
                -(tw3i * t11.im - tw5i * t10.im - tw2i * t9.im + tw1i * t8.im + tw4i * t7.im),
                tw3i * t11.re - tw5i * t10.re - tw2i * t9.re + tw1i * t8.re + tw4i * t7.re
            );
            (ch[0 + ido * (k + l1 * 3)], ch[0 + ido * (k + l1 * 8)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW4R * t2.re + TW3R * t3.re + TW1R * t4.re + TW5R * t5.re + TW2R * t6.re,
                t1.im + TW4R * t2.im + TW3R * t3.im + TW1R * t4.im + TW5R * t5.im + TW2R * t6.im
            );
            let cb = Complex::new(
                -(tw4i * t11.im - tw3i * t10.im + tw1i * t9.im + tw5i * t8.im - tw2i * t7.im),
                tw4i * t11.re - tw3i * t10.re + tw1i * t9.re + tw5i * t8.re - tw2i * t7.re
            );
            (ch[0 + ido * (k + l1 * 4)], ch[0 + ido * (k + l1 * 7)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW5R * t2.re + TW1R * t3.re + TW4R * t4.re + TW2R * t5.re + TW3R * t6.re,
                t1.im + TW5R * t2.im + TW1R * t3.im + TW4R * t4.im + TW2R * t5.im + TW3R * t6.im
            );
            let cb = Complex::new(
                -(tw5i * t11.im - tw1i * t10.im + tw4i * t9.im - tw2i * t8.im + tw3i * t7.im),
                tw5i * t11.re - tw1i * t10.re + tw4i * t9.re - tw2i * t8.re + tw3i * t7.re
            );
            (ch[0 + ido * (k + l1 * 5)], ch[0 + ido * (k + l1 * 6)]) = pmc(ca, cb);
        }
    }
    else {
        for k in 0..l1 {
            let t1 = cc[0 + ido * (0 + k * cdim)];
            let (t2, t11) = pmc(cc[0 + ido * (1 + k * cdim)], cc[0 + ido * (10 + k * cdim)]);
            let (t3, t10) = pmc(cc[0 + ido * (2 + k * cdim)], cc[0 + ido * (9 + k * cdim)]);
            let (t4, t9) = pmc(cc[0 + ido * (3 + k * cdim)], cc[0 + ido * (8 + k * cdim)]);
            let (t5, t8) = pmc(cc[0 + ido * (4 + k * cdim)], cc[0 + ido * (7 + k * cdim)]);
            let (t6, t7) = pmc(cc[0 + ido * (5 + k * cdim)], cc[0 + ido * (6 + k * cdim)]);

            ch[0 + ido * (k + l1 * 0)] = t1 + t2 + t3 + t4 + t5 + t6;

            let ca = Complex::new( 
                t1.re + TW1R * t2.re + TW2R * t3.re + TW3R * t4.re + TW4R * t5.re + TW5R * t6.re,
                t1.im + TW1R * t2.im + TW2R * t3.im + TW3R * t4.im + TW4R * t5.im + TW5R * t6.im
            );
            let cb = Complex::new(
                -(tw1i * t11.im + tw2i * t10.im + tw3i * t9.im + tw4i * t8.im + tw5i * t7.im),
                tw1i * t11.re + tw2i * t10.re + tw3i * t9.re + tw4i * t8.re + tw5i * t7.re
            );
            (ch[0 + ido * (k + l1 * 1)], ch[0 + ido * (k + l1 * 10)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW2R * t2.re + TW4R * t3.re + TW5R * t4.re + TW3R * t5.re + TW1R * t6.re,
                t1.im + TW2R * t2.im + TW4R * t3.im + TW5R * t4.im + TW3R * t5.im + TW1R * t6.im
            );
            let cb = Complex::new(
                -(tw2i * t11.im + tw4i * t10.im - tw5i * t9.im - tw3i * t8.im - tw1i * t7.im),
                tw2i * t11.re + tw4i * t10.re - tw5i * t9.re - tw3i * t8.re - tw1i * t7.re
            );
            (ch[0 + ido * (k + l1 * 2)], ch[0 + ido * (k + l1 * 9)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW3R * t2.re + TW5R * t3.re + TW2R * t4.re + TW1R * t5.re + TW4R * t6.re,
                t1.im + TW3R * t2.im + TW5R * t3.im + TW2R * t4.im + TW1R * t5.im + TW4R * t6.im
            );
            let cb = Complex::new(
                -(tw3i * t11.im - tw5i * t10.im - tw2i * t9.im + tw1i * t8.im + tw4i * t7.im),
                tw3i * t11.re - tw5i * t10.re - tw2i * t9.re + tw1i * t8.re + tw4i * t7.re
            );
            (ch[0 + ido * (k + l1 * 3)], ch[0 + ido * (k + l1 * 8)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW4R * t2.re + TW3R * t3.re + TW1R * t4.re + TW5R * t5.re + TW2R * t6.re,
                t1.im + TW4R * t2.im + TW3R * t3.im + TW1R * t4.im + TW5R * t5.im + TW2R * t6.im
            );
            let cb = Complex::new(
                -(tw4i * t11.im - tw3i * t10.im + tw1i * t9.im + tw5i * t8.im - tw2i * t7.im),
                tw4i * t11.re - tw3i * t10.re + tw1i * t9.re + tw5i * t8.re - tw2i * t7.re
            );
            (ch[0 + ido * (k + l1 * 4)], ch[0 + ido * (k + l1 * 7)]) = pmc(ca, cb);

            let ca = Complex::new(
                t1.re + TW5R * t2.re + TW1R * t3.re + TW4R * t4.re + TW2R * t5.re + TW3R * t6.re,
                t1.im + TW5R * t2.im + TW1R * t3.im + TW4R * t4.im + TW2R * t5.im + TW3R * t6.im
            );
            let cb = Complex::new(
                -(tw5i * t11.im - tw1i * t10.im + tw4i * t9.im - tw2i * t8.im + tw3i * t7.im),
                tw5i * t11.re - tw1i * t10.re + tw4i * t9.re - tw2i * t8.re + tw3i * t7.re
            );
            (ch[0 + ido * (k + l1 * 5)], ch[0 + ido * (k + l1 * 6)]) = pmc(ca, cb);

            for i in 1..ido {
                let t1 = cc[i + ido * (0 + k * cdim)];
                let (t2, t11) = pmc(cc[i + ido * (1 + k * cdim)], cc[i + ido * (10 + k * cdim)]);
                let (t3, t10) = pmc(cc[i + ido * (2 + k * cdim)], cc[i + ido * (9 + k * cdim)]);
                let (t4, t9) = pmc(cc[i + ido * (3 + k * cdim)], cc[i + ido * (8 + k * cdim)]);
                let (t5, t8) = pmc(cc[i + ido * (4 + k * cdim)], cc[i + ido * (7 + k * cdim)]);
                let (t6, t7) = pmc(cc[i + ido * (5 + k * cdim)], cc[i + ido * (6 + k * cdim)]);

                ch[i + ido * (k + l1 * 0)] = t1 + t2 + t3 + t4 + t5 + t6;

                let ca = Complex::new(
                    t1.re + TW1R * t2.re + TW2R * t3.re + TW3R * t4.re + TW4R * t5.re + TW5R * t6.re,
                    t1.im + TW1R * t2.im + TW2R * t3.im + TW3R * t4.im + TW4R * t5.im + TW5R * t6.im
                );
                let cb = Complex::new(
                    -(tw1i * t11.im + tw2i * t10.im + tw3i * t9.im + tw4i * t8.im + tw5i * t7.im),
                    tw1i * t11.re + tw2i * t10.re + tw3i * t9.re + tw4i * t8.re + tw5i * t7.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 1)] = if sign < 0 { wa[i - 1 + 0 * (ido - 1)].conj() * da } else { wa[i - 1 + 0 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 10)] = if sign < 0 { wa[i - 1 + 9 * (ido - 1)].conj() * db } else { wa[i - 1 + 9 * (ido - 1)] * db };

                let ca = Complex::new(
                    t1.re + TW2R * t2.re + TW4R * t3.re + TW5R * t4.re + TW3R * t5.re + TW1R * t6.re,
                    t1.im + TW2R * t2.im + TW4R * t3.im + TW5R * t4.im + TW3R * t5.im + TW1R * t6.im
                );
                let cb = Complex::new(
                    -(tw2i * t11.im + tw4i * t10.im - tw5i * t9.im - tw3i * t8.im - tw1i * t7.im),
                    tw2i * t11.re + tw4i * t10.re - tw5i * t9.re - tw3i * t8.re - tw1i * t7.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 2)] = if sign < 0 { wa[i - 1 + 1 * (ido - 1)].conj() * da } else { wa[i - 1 + 1 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 9)] = if sign < 0 { wa[i - 1 + 8 * (ido - 1)].conj() * db } else { wa[i - 1 + 8 * (ido - 1)] * db };

                let ca = Complex::new(
                    t1.re + TW3R * t2.re + TW5R * t3.re + TW2R * t4.re + TW1R * t5.re + TW4R * t6.re,
                    t1.im + TW3R * t2.im + TW5R * t3.im + TW2R * t4.im + TW1R * t5.im + TW4R * t6.im
                );
                let cb = Complex::new(
                    -(tw3i * t11.im - tw5i * t10.im - tw2i * t9.im + tw1i * t8.im + tw4i * t7.im),
                    tw3i * t11.re - tw5i * t10.re - tw2i * t9.re + tw1i * t8.re + tw4i * t7.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 3)] = if sign < 0 { wa[i - 1 + 2 * (ido - 1)].conj() * da } else { wa[i - 1 + 2 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 8)] = if sign < 0 { wa[i - 1 + 7 * (ido - 1)].conj() * db } else { wa[i - 1 + 7 * (ido - 1)] * db };

                let ca = Complex::new(
                    t1.re + TW4R * t2.re + TW3R * t3.re + TW1R * t4.re + TW5R * t5.re + TW2R * t6.re,
                    t1.im + TW4R * t2.im + TW3R * t3.im + TW1R * t4.im + TW5R * t5.im + TW2R * t6.im
                );
                let cb = Complex::new(
                    -(tw4i * t11.im - tw3i * t10.im + tw1i * t9.im + tw5i * t8.im - tw2i * t7.im),
                    tw4i * t11.re - tw3i * t10.re + tw1i * t9.re + tw5i * t8.re - tw2i * t7.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 4)] = if sign < 0 { wa[i - 1 + 3 * (ido - 1)].conj() * da } else { wa[i - 1 + 3 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 7)] = if sign < 0 { wa[i - 1 + 6 * (ido - 1)].conj() * db } else { wa[i - 1 + 6 * (ido - 1)] * db };

                let ca = Complex::new(
                    t1.re + TW5R * t2.re + TW1R * t3.re + TW4R * t4.re + TW2R * t5.re + TW3R * t6.re,
                    t1.im + TW5R * t2.im + TW1R * t3.im + TW4R * t4.im + TW2R * t5.im + TW3R * t6.im
                );
                let cb = Complex::new(
                    -(tw5i * t11.im - tw1i * t10.im + tw4i * t9.im - tw2i * t8.im + tw3i * t7.im),
                    tw5i * t11.re - tw1i * t10.re + tw4i * t9.re - tw2i * t8.re + tw3i * t7.re
                );
                let (da, db) = pmc(ca, cb);
                ch[i + ido * (k + l1 * 5)] = if sign < 0 { wa[i - 1 + 4 * (ido - 1)].conj() * da } else { wa[i - 1 + 4 * (ido - 1)] * da };
                ch[i + ido * (k + l1 * 6)] = if sign < 0 { wa[i - 1 + 5 * (ido - 1)].conj() * db } else { wa[i - 1 + 5 * (ido - 1)] * db };
            }
        }
    }
}

fn passg(ido: usize, ip: usize, l1: usize, cc: &mut [Complex], ch: &mut [Complex], wa: &[Complex], csarr: &[Complex], sign: i32) {
    let cdim = ip;
    let ipph = (ip + 1) >> 1;
    let idl1 = ido * l1;

    let mut wal = vec![Complex::new(0.0, 0.0); ip];
    wal[0] = Complex::new(1.0, 0.0);
    for i in 1..ip { wal[i] = Complex::new(csarr[i].re, sign as f64 * csarr[i].im); }

    for k in 0..l1 {
        for i in 0..ido {
            ch[i + ido * (k + l1 * 0)] = cc[i + ido * (0 + k * cdim)];
        }
    }

    for j in 1..ipph {
        let jc = ip - j;
        for k in 0..l1 {
            for i in 0..ido {
                (ch[i + ido * (k + l1 * j)], ch[i + ido * (k + l1 * jc)]) = pmc(cc[i + ido * (j + k * cdim)], cc[i + ido * (jc + k * cdim)]);
            }
        }
    }

    for k in 0..l1 {
        for i in 0..ido {
            let mut sum = ch[i + ido * (k + l1 * 0)];
            for j in 1..ipph {
                sum += ch[i + ido * (k + l1 * j)];
            }
            cc[i + ido * (k + l1 * 0)] = sum;
        }
    }

    for l in 1..ipph {
        let lc = ip - l;

        for ik in 0..idl1 {
            cc[ik + idl1 * l] = Complex::new(
                ch[ik + idl1 * 0].re + wal[l].re * ch[ik + idl1 * 1].re + wal[2 * l].re * ch[ik + idl1 * 2].re,
                ch[ik + idl1 * 0].im + wal[l].re * ch[ik + idl1 * 1].im + wal[2 * l].re * ch[ik + idl1 * 2].im
            );
            cc[ik + idl1 * lc] = Complex::new(
                -wal[l].im * ch[ik + idl1 * (ip - 1)].im - wal[2 * l].im * ch[ik + idl1 * (ip - 2)].im,
                wal[l].im * ch[ik + idl1 * (ip - 1)].re + wal[2 * l].im * ch[ik + idl1 * (ip - 2)].re
            );
        }

        let mut iwal = 2 * l;
        let mut j = 3;
        let mut jc = ip - 3;

        while j < ipph - 1 {
            iwal += l;
            if iwal > ip { iwal -= ip; }
            let xwal = wal[iwal];

            iwal += l;
            if iwal > ip { iwal -= ip; }
            let xwal2 = wal[iwal];

            for ik in 0..idl1 {
                cc[ik + idl1 * l].re += ch[ik + idl1 * j].re * xwal.re + ch[ik + idl1 * (j + 1)].re * xwal2.re;
                cc[ik + idl1 * l].im += ch[ik + idl1 * j].im * xwal.re + ch[ik + idl1 * (j + 1)].im * xwal2.re;
                cc[ik + idl1 * lc].re -= ch[ik + idl1 * jc].im * xwal.im + ch[ik + idl1 * (jc - 1)].im * xwal2.im;
                cc[ik + idl1 * lc].im += ch[ik + idl1 * jc].re * xwal.im + ch[ik + idl1 * (jc - 1)].re * xwal2.im;
            }
            j += 2;
            jc -= 2;
        }

        while j < ipph {
            iwal += l;
            if iwal > ip { iwal -= ip; }
            let xwal = wal[iwal];

            for ik in 0..idl1 {
                cc[ik + idl1 * l].re += ch[ik + idl1 * j].re * xwal.re;
                cc[ik + idl1 * l].im += ch[ik + idl1 * j].im * xwal.re;
                cc[ik + idl1 * lc].re -= ch[ik + idl1 * jc].im * xwal.im;
                cc[ik + idl1 * lc].im += ch[ik + idl1 * jc].re * xwal.im;
            }
            j += 1;
            jc -= 1;
        }
    }

    if ido == 1 {
        for j in 1..ipph {
        let jc = ip - j;
            for ik in 0..idl1 {
                let t1 = cc[ik + idl1 * j];
                let t2 = cc[ik + idl1 * jc];
                (cc[ik + idl1 * j], cc[ik + idl1 * jc]) = pmc(t1, t2);
            }
        }
    }
    else {
        for j in 1..ipph {
            let jc = ip - j;
            for k in 0..l1 {
                let t1 = cc[0 + ido * (k + l1 * j)];
                let t2 = cc[0 + ido * (k + l1 * jc)];
                (cc[0 + ido * (k + l1 * j)], cc[0 + ido * (k + l1 * jc)]) = pmc(t1, t2);

                for i in 1..ido {
                    let (x1, x2) = pmc(cc[i + ido * (k + l1 * j)], cc[i + ido * (k + l1 * jc)]);
                    let idij = (j - 1) * (ido - 1) + i - 1;
                    let idij2 = (jc - 1) * (ido - 1) + i - 1;
                    cc[i + ido * (k + l1 * j)] = if sign < 0 { wa[idij].conj() * x1 } else { wa[idij] * x1 };
                    cc[i + ido * (k + l1 * jc)] = if sign < 0 { wa[idij2].conj() * x2 } else { wa[idij2] * x2 };
                }
            }
        }
    }
}

#[derive(Clone)]
pub struct FactorData {
    pub fct: usize,
    pub tw: Vec<Complex>,
    pub tws: Vec<Complex>,
}

pub struct CooleyTukey {
    len: usize,
    fct: Vec<FactorData>,
}

impl CooleyTukey {
    pub fn new(len: usize) -> CooleyTukey {
        let mut plan = CooleyTukey { len, fct: Vec::new() };
        if len < 2 { return plan; }

        plan.factorize();
        plan.compute_twiddle();
        return plan;
    }

    fn factorize(&mut self) {
        let mut len = self.len;

        while len & 3 == 0 {
            self.fct.push(FactorData { fct: 4, tw: Vec::new(), tws: Vec::new() });
            len >>= 2;
        }

        if len & 1 == 0 {
            len >>= 1;
            self.fct.push(FactorData { fct: 2, tw: Vec::new(), tws: Vec::new() });
            let fctlen = self.fct.len(); self.fct.swap(0, fctlen - 1);
        }

        let mut maxl = (len as f64).sqrt() as usize + 1;
        let mut divisor = 3;
        while len > 1 && divisor < maxl {
            if len % divisor == 0 {
                while len % divisor == 0 {
                    self.fct.push(FactorData { fct: divisor, tw: Vec::new(), tws: Vec::new() });
                    len /= divisor;
                }
                maxl = (len as f64).sqrt() as usize + 1;
            }
            divisor += 2;
        }

        if len > 1 { self.fct.push(FactorData { fct: len, tw: Vec::new(), tws: Vec::new() }); }
    }

    fn compute_twiddle(&mut self) {
        let len = self.len;
        let mut twid = vec![Complex::new(0.0, 0.0); len];
        sincos_2pibyn(len, &mut twid);

        let mut l1 = 1;

        for k in 0..self.fct.len() {
            let ip = self.fct[k].fct;
            let ido = len / (l1 * ip);

            let tw_size = (ip - 1) * (ido - 1);
            self.fct[k].tw = vec![Complex::new(0.0, 0.0); tw_size];

            for j in 1..ip {
                for i in 1..ido {
                    let idx = (j - 1) * (ido - 1) + i - 1;
                    self.fct[k].tw[idx] = twid[j * l1 * i];
                }
            }

            if ip > 11 {
                self.fct[k].tws = vec![Complex::new(0.0, 0.0); ip];
                for j in 0..ip { self.fct[k].tws[j] = twid[j * l1 * ido]; }
            }

            l1 *= ip;
        }
    }

    pub fn forward(&self, data: &mut [Complex], fct: f64) { self.fft(data, fct, -1); }
    pub fn backward(&self, data: &mut [Complex], fct: f64) { self.fft(data, fct, 1); }

    fn fft(&self, data: &mut [Complex], fct: f64, sign: i32) {
        if self.len != data.len() { panic!("Provided buffer({}) not compatible with CooleyTukey({}).", data.len(), self.len); }
        if self.len < 2 { return; }

        let mut l1 = 1;
        let mut ch = vec![Complex::new(0.0, 0.0); data.len()];
        let (mut p1, mut p2) = (&mut data[..], &mut ch[..]);

        for k1 in 0..self.fct.len() {
            let ip = self.fct[k1].fct;
            let l2 = ip * l1;
            let ido = self.len / l2;

            match ip {
                4 => pass4(ido, l1, p1, p2, &self.fct[k1].tw, sign),
                2 => pass2(ido, l1, p1, p2, &self.fct[k1].tw, sign),
                3 => pass3(ido, l1, p1, p2, &self.fct[k1].tw, sign),
                5 => pass5(ido, l1, p1, p2, &self.fct[k1].tw, sign),
                7 => pass7(ido, l1, p1, p2, &self.fct[k1].tw, sign),
                11 => pass11(ido, l1, p1, p2, &self.fct[k1].tw, sign),
                _ => { passg(ido, ip, l1, p1, p2, &self.fct[k1].tw, &self.fct[k1].tws, sign); (p1, p2) = (p2, p1); }
            }

            (p1, p2, l1) = (p2, p1, l2);
        }

        if fct != 1.0 { p1.iter_mut().for_each(|d| *d *= fct); }
        if p1.as_ptr() != data.as_ptr() { data.copy_from_slice(&ch); }
    }

    pub fn len(&self) -> usize { self.len }
}