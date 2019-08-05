#![allow(dead_code)]

mod cov;
use crate::cov::*;

#[derive(Debug)]
struct XMeas (Vec3, Cov3);
// struct XMeas (f64,  f64);

impl XMeas {
    fn to_string (&self) -> String {
        // let XMeas(v, cv) = self;
        String::from("XMeas!!")
    }

}

fn main() {
    let vxc3: [f64; 6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let xc3: Cov3 = Cov3 {v: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]};
    let x3: XMeas = XMeas (
        Vec3 ( [1.0, 2.0, 3.0] ),
        xc3,
    );
    println!("vxc3 is {:?}", vxc3);
    println!("x3.1 is {}", x3.1.to_string());
    println!("x3 is {}", x3.to_string());
    println!("x3 is also {:?}", x3);
    println!("x3 is also {:#?}", x3);
}
