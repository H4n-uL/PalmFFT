pub mod bluestein; pub mod cooleytukey;
use self::{bluestein::Bluestein, cooleytukey::CooleyTukey};
use crate::{
    math::{cost_guess, good_size, largest_prime_factor},
    Complex, Result
};

pub enum CfftPlan {
    Ct(CooleyTukey),
    Bs(Bluestein)
}

impl CfftPlan {
    pub fn new(length: usize) -> Self {
        if length < 50 || largest_prime_factor(length) <= libm::sqrt(length as f64) as usize {
            return Self::Ct(CooleyTukey::new(length));
        }
        let ct_cost = cost_guess(length);
        let bs_cost = 3.0 * cost_guess(good_size(2 * length - 1));

        if bs_cost < ct_cost {
            return Self::Bs(Bluestein::new(length));
        }
        return Self::Ct(CooleyTukey::new(length));
    }

    pub fn forward(&self, data: &mut [Complex], fct: f64) -> Result {
        match self {
            Self::Ct(ct) => { ct.forward(data, fct) }
            Self::Bs(bs) => { bs.forward(data, fct) }
        }
    }

    pub fn backward(&self, data: &mut [Complex], fct: f64) -> Result {
        match self {
            Self::Ct(ct) => { ct.backward(data, fct) }
            Self::Bs(bs) => { bs.backward(data, fct) }
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Self::Ct(ct) => { ct.len() }
            Self::Bs(bs) => { bs.len() }
        }
    }
}