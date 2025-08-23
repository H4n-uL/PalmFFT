mod complex; mod algorithms; mod math;
pub use complex::Complex;
pub use algorithms::CfftPlan;

#[cfg(feature = "stddep")]
pub use num_complex::Complex64 as NumComplex;

#[cfg(feature = "stddep")]
impl From<NumComplex> for Complex {
    #[inline] fn from(c: NumComplex) -> Self { Self::new(c.re, c.im) }
}

#[cfg(feature = "stddep")]
impl From<Complex> for NumComplex {
    #[inline] fn from(c: Complex) -> Self { Self::new(c.re, c.im) }
}