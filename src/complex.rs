use core::ops::{Neg, Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};
use std::fmt::{Debug, Display};

#[derive(Clone, Copy, PartialEq)]
pub struct Complex { pub re: f64, pub im: f64 }

impl Complex {
    #[inline] pub fn new(re: f64, im: f64) -> Self { Self { re, im } }
    #[inline] pub fn zero() -> Self { Self { re: 0.0, im: 0.0 } }
    #[inline] pub fn one() -> Self { Self { re: 1.0, im: 0.0 } }
    #[inline] pub fn i() -> Self { Self { re: 0.0, im: 1.0 } }

    #[inline] pub fn conj(&self) -> Self { Self { re: self.re, im: -self.im } } // conjugate
    #[inline] pub fn rot90(&self) -> Self { Self { re: -self.im, im: self.re } } // rotate 90 degrees
    #[inline] pub fn rotm90(&self) -> Self { Self { re: self.im, im: -self.re } } // rotate -90 degrees
    #[inline] pub fn dot(&self, rhs: Complex) -> f64 { self.re * rhs.re + self.im * rhs.im } // dot product
    #[inline] pub fn norm_sqr(&self) -> f64 { self.dot(*self) } // complex norm squared
    #[inline] pub fn norm(&self) -> f64 { self.norm_sqr().sqrt() } // complex norm
    #[inline] pub fn arg(&self) -> f64 { self.im.atan2(self.re) } // complex argument

    #[inline] pub fn to_polar(&self) -> (f64, f64) { (self.norm(), self.arg()) } // convert to polar form
    #[inline] pub fn from_polar(re: f64, theta: f64) -> Self { Self { re: theta.cos(), im: theta.sin() } * re } // convert from polar form

    #[inline] pub fn exp(&self) -> Self { Complex::from_polar(self.re.exp(), self.im) } // complex natural exponentation
    #[inline] pub fn ln(&self) -> Self { Self { re: self.norm().ln(), im: self.arg() } } // complex natural logarithm

    #[inline] pub fn powf(&self, rhs: f64) -> Self { Complex::from_polar(self.norm().powf(rhs), rhs * self.arg()) }
    #[inline] pub fn powc(&self, rhs: Complex) -> Self { Complex::from_polar(rhs.dot(self.ln().conj()).exp(), rhs.dot(self.ln().conj().rot90())) } 
    #[inline] pub fn rootf(&self, rhs: f64) -> Self { self.powf(1.0 / rhs) }
    #[inline] pub fn rootc(&self, rhs: Complex) -> Self { self.powc(1.0 / rhs) }

    #[inline] pub fn is_nan(self) -> bool { self.re.is_nan() || self.im.is_nan() }
    #[inline] pub fn is_infinite(self) -> bool { !self.is_nan() && (self.re.is_infinite() || self.im.is_infinite()) }
    #[inline] pub fn is_finite(self) -> bool { self.re.is_finite() && self.im.is_finite() }
    #[inline] pub fn is_normal(self) -> bool { self.re.is_normal() && self.im.is_normal() }
    #[inline] pub fn is_subnormal(self) -> bool { self.re.is_subnormal() || self.im.is_subnormal() }
}

impl Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (self.re, self.im) {
            (0.0, 0.0) => write!(f, "0"),
            (0.0, _) => write!(f, "{}i", self.im),
            (_, 0.0) => write!(f, "{}", self.re),
            (_, _) => write!(f, "{}{}{}i", self.re, if self.im < 0.0 { "" } else { "+" }, self.im)
        }
    }
}

impl Debug for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (self.re, self.im) {
            (0.0, 0.0) => write!(f, "0"),
            (0.0, _) => write!(f, "{}i", self.im),
            (_, 0.0) => write!(f, "{}", self.re),
            (_, _) => write!(f, "{}{}{}i", self.re, if self.im < 0.0 { "" } else { "+" }, self.im)
        }
    }
}

impl Neg for Complex { type Output = Self; #[inline] fn neg(self) -> Self { Self { re: -self.re, im: -self.im } } }

impl Add for Complex { type Output = Self; #[inline] fn add(self, rhs: Complex) -> Self { Self { re: self.re + rhs.re, im: self.im + rhs.im } } }
impl Add<f64> for Complex { type Output = Self; #[inline] fn add(self, rhs: f64) -> Self { Self { re: self.re + rhs, im: self.im } } }
impl Add<Complex> for f64 { type Output = Complex; #[inline] fn add(self, rhs: Complex) -> Complex { rhs + self } }
impl AddAssign for Complex { #[inline] fn add_assign(&mut self, rhs: Complex) { *self = *self + rhs; } }
impl AddAssign<f64> for Complex { #[inline] fn add_assign(&mut self, rhs: f64) { *self = *self + rhs; } }

impl Sub for Complex { type Output = Self; #[inline] fn sub(self, rhs: Complex) -> Self { Self { re: self.re - rhs.re, im: self.im - rhs.im } } }
impl Sub<f64> for Complex { type Output = Self; #[inline] fn sub(self, rhs: f64) -> Self { Self { re: self.re - rhs, im: self.im } } }
impl Sub<Complex> for f64 { type Output = Complex; #[inline] fn sub(self, rhs: Complex) -> Complex { rhs - self } }
impl SubAssign for Complex { #[inline] fn sub_assign(&mut self, rhs: Complex) { *self = *self - rhs; } }
impl SubAssign<f64> for Complex { #[inline] fn sub_assign(&mut self, rhs: f64) { *self = *self - rhs; } }

impl Mul for Complex { type Output = Self; #[inline] fn mul(self, rhs: Complex) -> Self { Self { re: self.dot(rhs.conj()), im: self.dot(rhs.rotm90().conj()) } } }
impl Mul<f64> for Complex { type Output = Self; #[inline] fn mul(self, rhs: f64) -> Self { Self { re: self.re * rhs, im: self.im * rhs } } }
impl Mul<Complex> for f64 { type Output = Complex; #[inline] fn mul(self, rhs: Complex) -> Complex { rhs * self } }
impl MulAssign for Complex { #[inline] fn mul_assign(&mut self, rhs: Complex) { *self = *self * rhs; } }
impl MulAssign<f64> for Complex { #[inline] fn mul_assign(&mut self, rhs: f64) { *self = *self * rhs; } }

impl Div for Complex { type Output = Self; #[inline] fn div(self, rhs: Complex) -> Self { Self { re: self.dot(rhs), im: self.dot(rhs.rot90()) } / rhs.norm_sqr() } }
impl Div<f64> for Complex { type Output = Self; #[inline] fn div(self, rhs: f64) -> Self { Self { re: self.re / rhs, im: self.im / rhs } } }
impl Div<Complex> for f64 { type Output = Complex; #[inline] fn div(self, rhs: Complex) -> Complex { self * rhs.conj() / rhs.norm() } }
impl DivAssign for Complex { #[inline] fn div_assign(&mut self, rhs: Complex) { *self = *self / rhs; } }
impl DivAssign<f64> for Complex { #[inline] fn div_assign(&mut self, rhs: f64) { *self = *self / rhs; } }

impl From<f64> for Complex { #[inline] fn from(re: f64) -> Self { Self { re, im: 0.0 } } }