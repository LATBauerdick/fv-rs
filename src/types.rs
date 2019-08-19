
use crate::cov::*;

use std::fmt;

pub type Number = f64;

#[derive(Debug)]
pub struct XMeas(pub Vec3, pub Cov3);
impl XMeas {

}
impl fmt::Display for XMeas {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let XMeas(v, cv) = self;
        write!(f, "XMeas:{}\n{}", v, cv)
    }
}

#[derive(Debug)]
pub struct HMeas(pub Vec5, pub Cov5, pub Number);
impl HMeas {

}

#[derive(Debug)]
pub struct VHMeas {
    pub vertex:  XMeas,
    pub helices: Vec<HMeas>,
}

