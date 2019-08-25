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



// showMomentum :: HMeas -> String
// showMomentum h = "pt,pz,fi,E ->" <> (show <<< fromHMeas) h

// fromHMeas :: HMeas -> QMeas -- just drop the d0, z0 part... fix!!!!
// fromHMeas (HMeas h ch w2pt) = QMeas q cq w2pt where
//   q = subm 3 h
//   cq = subm2 3 ch

// -----------------------------------------------
// -- QMeas
// -- 3-vector and covariance matrix for momentum measurement
// --
// use std::f64::*;
// data QMeas = QMeas Vec3 Cov3 Number
// instance Show QMeas where
//   show = showQMeas
// -- print QMeas as a 4-momentum vector with errors, use pt and pz
// showQMeas :: QMeas -> String
use crate::inp::h_slurp;
#[test]
fn test_fvt() {
    let ds = std::fs::read_to_string("dat/tr05129e001412.dat").unwrap();
    let VHMeas {vertex: x, helices: hel} = h_slurp(ds).unwrap();
    for h in &hel { println!("{}", h) };
    for h in hel { println!("{}", QMeas::from(&h)) };
    println!("initial vertex position -> {}", x);
    let res = String::from("all good?");
    assert!( 1 == 1, "test failed with '{}'", res);

}
// testFVT l5 vm = do
//   let hel = helices vm
//   traverse_ (putStrLn <<< showHelix) hel
//   traverse_ (putStrLn <<< showMomentum) hel
//   doFitTest vm l5
//   putStrLn $ showProng <<< fit <<< hFilter l5 <<< vBlowup 10000.0 $ vm
//   pure ()

