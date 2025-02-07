use core::ops::{Neg, Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};
use std::fmt::{Debug, Display};

#[derive(Clone, Copy, PartialEq)]
pub struct Complex { pub r: f64, pub i: f64 }

impl Complex {
    pub fn new(r: f64, i: f64) -> Complex { Complex { r, i } }
    pub fn conj(&self) -> Complex { Complex { r: self.r, i: -self.i } } // conjugate
    pub fn rot90(&self) -> Complex { Complex { r: -self.i, i: self.r } } // rotate 90 degrees
    pub fn rotm90(&self) -> Complex { Complex { r: self.i, i: -self.r } } // rotate -90 degrees
    pub fn dot(&self, rhs: Complex) -> f64 { self.r * rhs.r + self.i * rhs.i } // dot product
    pub fn norm_sqr(&self) -> f64 { self.dot(*self) } // complex norm squared
    pub fn norm(&self) -> f64 { self.norm_sqr().sqrt() } // complex norm
    pub fn arg(&self) -> f64 { self.i.atan2(self.r) } // complex argument

    pub fn to_polar(&self) -> (f64, f64) { (self.norm(), self.arg()) } // convert to polar form
    pub fn from_polar(r: f64, theta: f64) -> Complex { Complex { r: theta.cos(), i: theta.sin() } * r } // convert from polar form

    pub fn exp(&self) -> Complex { Complex::from_polar(self.r.exp(), self.i) } // complex natural exponentation
    pub fn ln(&self) -> Complex { Complex { r: self.norm().ln(), i: self.arg() } } // complex natural logarithm

    pub fn powf(&self, rhs: f64) -> Complex { Complex::from_polar(self.norm().powf(rhs), rhs * self.arg()) }
    pub fn powc(&self, rhs: Complex) -> Complex { Complex::from_polar(rhs.dot(self.ln().conj()).exp(), rhs.dot(self.ln().conj().rot90())) } 
    pub fn rootf(&self, rhs: f64) -> Complex { self.powf(1.0 / rhs) }
    pub fn rootc(&self, rhs: Complex) -> Complex { self.powc(1.0 / rhs) }

    pub fn is_nan(self) -> bool { self.r.is_nan() || self.i.is_nan() }
    pub fn is_infinite(self) -> bool { !self.is_nan() && (self.r.is_infinite() || self.i.is_infinite()) }
    pub fn is_finite(self) -> bool { self.r.is_finite() && self.i.is_finite() }
    pub fn is_normal(self) -> bool { self.r.is_normal() && self.i.is_normal() }
    pub fn is_subnormal(self) -> bool { self.r.is_subnormal() || self.i.is_subnormal() }
}

impl Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (self.r, self.i) {
            (0.0, 0.0) => write!(f, "0"),
            (0.0, _) => write!(f, "{}i", self.i),
            (_, 0.0) => write!(f, "{}", self.r),
            (_, _) => write!(f, "{}{}{}i", self.r, if self.i < 0.0 { "" } else { "+" }, self.i)
        }
    }
}

impl Debug for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (self.r, self.i) {
            (0.0, 0.0) => write!(f, "0"),
            (0.0, _) => write!(f, "{}i", self.i),
            (_, 0.0) => write!(f, "{}", self.r),
            (_, _) => write!(f, "{}{}{}i", self.r, if self.i < 0.0 { "" } else { "+" }, self.i)
        }
    }
}

impl Neg for Complex { type Output = Complex; fn neg(self) -> Complex { Complex { r: -self.r, i: -self.i } } }

impl Add for Complex { type Output = Complex; fn add(self, rhs: Complex) -> Complex { Complex { r: self.r + rhs.r, i: self.i + rhs.i } } }
impl Add<f64> for Complex { type Output = Complex; fn add(self, rhs: f64) -> Complex { Complex { r: self.r + rhs, i: self.i } } }
impl Add<Complex> for f64 { type Output = Complex; fn add(self, rhs: Complex) -> Complex { rhs + self } }
impl AddAssign for Complex { fn add_assign(&mut self, rhs: Complex) { *self = *self + rhs; } }
impl AddAssign<f64> for Complex { fn add_assign(&mut self, rhs: f64) { *self = *self + rhs; } }

impl Sub for Complex { type Output = Complex; fn sub(self, rhs: Complex) -> Complex { Complex { r: self.r - rhs.r, i: self.i - rhs.i } } }
impl Sub<f64> for Complex { type Output = Complex; fn sub(self, rhs: f64) -> Complex { Complex { r: self.r - rhs, i: self.i } } }
impl Sub<Complex> for f64 { type Output = Complex; fn sub(self, rhs: Complex) -> Complex { rhs - self } }
impl SubAssign for Complex { fn sub_assign(&mut self, rhs: Complex) { *self = *self - rhs; } }
impl SubAssign<f64> for Complex { fn sub_assign(&mut self, rhs: f64) { *self = *self - rhs; } }

impl Mul for Complex { type Output = Complex; fn mul(self, rhs: Complex) -> Complex { Complex { r: self.dot(rhs.conj()), i: self.dot(rhs.rotm90().conj()) } } }
impl Mul<f64> for Complex { type Output = Complex; fn mul(self, rhs: f64) -> Complex { Complex { r: self.r * rhs, i: self.i * rhs } } }
impl Mul<Complex> for f64 { type Output = Complex; fn mul(self, rhs: Complex) -> Complex { rhs * self } }
impl MulAssign for Complex { fn mul_assign(&mut self, rhs: Complex) { *self = *self * rhs; } }
impl MulAssign<f64> for Complex { fn mul_assign(&mut self, rhs: f64) { *self = *self * rhs; } }

impl Div for Complex { type Output = Complex; fn div(self, rhs: Complex) -> Complex { Complex { r: self.dot(rhs), i: self.dot(rhs.rot90()) } / rhs.norm_sqr() } }
impl Div<f64> for Complex { type Output = Complex; fn div(self, rhs: f64) -> Complex { Complex { r: self.r / rhs, i: self.i / rhs } } }
impl Div<Complex> for f64 { type Output = Complex; fn div(self, rhs: Complex) -> Complex { self * rhs.conj() / rhs.norm() } }
impl DivAssign for Complex { fn div_assign(&mut self, rhs: Complex) { *self = *self / rhs; } }
impl DivAssign<f64> for Complex { fn div_assign(&mut self, rhs: f64) { *self = *self / rhs; } }