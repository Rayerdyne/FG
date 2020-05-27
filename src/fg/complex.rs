use std::fmt;

/// Holds a complex number
#[derive(Debug, Clone, Copy)]
pub struct Complex {
    /// The real part
    pub re: f64,
    /// The imaginary part
    pub im: f64,
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + j * {})", self.re, self.im)?;
        Ok(())
    }
}

impl std::ops::Add for Complex {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {  re: self.re + other.re,
                im: self.im + other.im  }
    }
}
impl std::ops::AddAssign for Complex {
    fn add_assign(&mut self, c: Self) {  self.re += c.re;
                                         self.im += c.im;    }
}
impl std::ops::Sub for Complex {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {  re: self.re - other.re,
                im: self.im - other.im  }
    }
}
impl std::ops::SubAssign for Complex {
    fn sub_assign(&mut self, c: Self) { self.re -= c.re;
                                        self.im -= c.im;    }
}
impl std::ops::Mul for Complex {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self  {     re: self.re * other.re - self.im * other.im,
                    im: self.re * other.im + self.im * other.re }
    }
}
impl std::ops::MulAssign for Complex {
    fn mul_assign(&mut self, other: Self) {
        let new_re = self.re * other.re - self.im * other.im;
        let new_im = self.re * other.im + self.im * other.re;
        self.re = new_re;
        self.im = new_im;
    }
}
impl std::ops::MulAssign<f64> for Complex {
    fn mul_assign(&mut self, div: f64) {    self.re *= div;
                                            self.im *= div;     }
}
impl std::ops::DivAssign<f64> for Complex {
    fn div_assign(&mut self, div: f64) {    self.re /= div;
                                            self.im /= div;     }
}

impl Complex {
    /// Multiplies the complex by the imaginary unit.
    pub fn times_j(&mut self) -> Self {
        let temp = self.re;
        self.re = -self.im;
        self.im = temp;
        *self
    }
    /// Returns a complex containing 0 + 0 * j
    pub fn zero() -> Self {
        Complex {   re: 0.0,
                    im: 0.0     }
    }
    /// Multiplies the complex by e^{j * theta}
    pub fn expj(theta: f64) -> Self {
        Complex {   re: theta.cos(),
                    im: theta.sin()     }
    }
}