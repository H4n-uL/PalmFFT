use crate::{
    algorithms::cooleytukey::CooleyTukey,
    complex::Complex,
    math::{good_size, sincos_2pibyn}
};

pub struct Bluestein {
    n: usize,
    n2: usize,
    plan: CooleyTukey,
    bk: Vec<Complex>,
    bkf: Vec<Complex>,
}

impl Bluestein {
    pub fn new(length: usize) -> Bluestein {
        let n = length;
        let n2 = if length == 0 { 0 } else { good_size(n * 2 - 1) };
        let mut plan = Bluestein {
            n, n2, plan: CooleyTukey::new(n2),
            bk: vec![Complex::new(0.0, 0.0); n],
            bkf: vec![Complex::new(0.0, 0.0); n2]
        };

        if plan.n < 2 { return plan; }

        let mut tmp = vec![Complex::new(0.0, 0.0); n * 2];
        sincos_2pibyn(n * 2, &mut tmp);
        plan.bk[0] = Complex::new(1.0, 0.0);

        let mut coeff = 0;
        for m in 1..n {
            coeff += 2 * m - 1;
            if coeff >= 2 * n { coeff -= 2 * n; }
            plan.bk[m] = tmp[coeff];
        }

        let xn2 = 1.0 / (n2 as f64);
        plan.bkf[0] = plan.bk[0] * xn2;

        for m in 1..n {
            let norm = plan.bk[m] * xn2;
            (plan.bkf[m], plan.bkf[n2 - m]) = (norm, norm);
        }

        for m in n..=(n2 - n) { plan.bkf[m] = Complex::new(0.0, 0.0); }
        plan.plan.forward(&mut plan.bkf, 1.0);
        return plan;
    }

    pub fn forward(&self, data: &mut [Complex], fct: f64) { self.fft(data, fct, -1); }
    pub fn backward(&self, data: &mut [Complex], fct: f64) { self.fft(data, fct, 1); }

    pub fn fft(&self, data: &mut [Complex], fct: f64, sign: i32) {
        if self.n != data.len() { panic!("Provided buffer({}) not compatible with Bluestein({}).", data.len(), self.n); }
        if self.n < 2 { return; }
        let mut akf = vec![Complex::new(0.0, 0.0); self.n2];

        for m in 0..self.n { akf[m] = data[m] * if sign > 0 { self.bk[m] } else { self.bk[m].conj() }; }
        for m in self.n..self.n2 { akf[m] = Complex::new(0.0, 0.0); }

        self.plan.forward(&mut akf, fct);
        for m in 0..self.n2 { akf[m] *= if sign > 0 { self.bkf[m].conj() } else { self.bkf[m] }; }
        self.plan.backward(&mut akf, fct);
        for m in 0..self.n { data[m] = akf[m] * if sign > 0 { self.bk[m] } else { self.bk[m].conj() }; }
    }

    pub fn len(&self) -> usize { self.n }
}