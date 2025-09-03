pub mod bluestein; pub mod cooleytukey;
use self::{bluestein::Bluestein, cooleytukey::CooleyTukey};
use crate::{
    math::{cost_guess, good_size, largest_prime_factor},
    Complex
};

pub struct CfftPlan {
    ct: Option<CooleyTukey>,
    bs: Option<Bluestein>
}

impl CfftPlan {
    pub fn new(length: usize) -> CfftPlan {
        let mut plan = CfftPlan { bs: None, ct: None };
        if length < 50 || largest_prime_factor(length) <= (length as f64).sqrt() as usize {
            plan.ct = Some(CooleyTukey::new(length));
            return plan;
        }
        let ct_cost = cost_guess(length);
        let bs_cost = 3.0 * cost_guess(good_size(2 * length - 1));

        if bs_cost < ct_cost { plan.bs = Some(Bluestein::new(length)); }
        else { plan.ct = Some(CooleyTukey::new(length)); }

        return plan;
    }

    pub fn forward(&self, data: &mut [Complex], fct: f64) {
        if let Some(ref ct) = self.ct { ct.forward(data, fct); }
        else if let Some(ref bs) = self.bs { bs.forward(data, fct); }
    }

    pub fn backward(&self, data: &mut [Complex], fct: f64) {
        if let Some(ref ct) = self.ct { ct.backward(data, fct); }
        else if let Some(ref bs) = self.bs { bs.backward(data, fct); }
    }

    pub fn len(&self) -> usize {
        if let Some(ref ct) = self.ct { ct.len() }
        else if let Some(ref bs) = self.bs { bs.len() }
        else { 0 }
    }
}