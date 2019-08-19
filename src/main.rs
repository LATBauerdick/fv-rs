#![allow(dead_code)]

use crate::types::*;

mod types;
mod cov;
mod chol;
mod inp;

//use crate::types::*;
use crate::cov::*;


fn main() {
    let vxc3: [Number; 6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let xc3 = Cov3::from([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let vc3: Vec3 = Vec3 { v: [1.0, 2.0, 3.0] };
    let x3: XMeas = XMeas(vc3, xc3);
    println!("vxc3 is {:?}", vxc3);
    println!("x3.1 is {}", x3.1);
    println!("x3 is {}", x3);
    println!("x3 is also {:?}", x3);
    println!("x3 is also {:#?}", x3);
}
