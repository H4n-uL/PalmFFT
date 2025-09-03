mod algorithms; mod math;
pub use algorithms::CfftPlan;

#[cfg(feature = "num-complex")]
pub use num_complex::Complex64 as Complex;

#[cfg(not(feature = "num-complex"))]
mod complex;
#[cfg(not(feature = "num-complex"))]
pub use complex::Complex;

pub trait ComplexExt {
    fn rot90(&self) -> Self;
    fn rotm90(&self) -> Self;
}

impl ComplexExt for Complex {
    #[inline]
    fn rot90(&self) -> Self {
        Complex::new(-self.im, self.re)
    }
    
    #[inline]
    fn rotm90(&self) -> Self {
        Complex::new(self.im, -self.re)
    }
}