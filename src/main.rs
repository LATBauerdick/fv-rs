#![allow(dead_code)]

use crate::types::*;

mod types;
mod fit;
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

#[test]
fn test_fvt() {
use crate::inp::h_slurp;
use crate::fit::*;
    println!("test_fvt-------------------------------------------------");
    let ds = std::fs::read_to_string("dat/tr05129e001412.dat").unwrap();
    let VHMeas {vertex: x, helices: hel} = h_slurp(ds).unwrap();
//   doFitTest vm l5
    let vm = VHMeas {vertex: x.blowup(10000.0), helices: hel};
    let l5 = vec![0_usize,2,3,4,5];

    for h in &vm.helices { println!("{}", h) };
    for h in &vm.helices { println!("{}", QMeas::from(h)) };

    println!("initial vertex position -> {}", vm.vertex);

    // for h in &vm.helices { println!("{}", PMeas::from(&QMeas::from(h))) };

    let mm = inv_mass( &vm.helices
                        .iter()
                        .map( |h| PMeas::from(&QMeas::from(h)) )
                        .collect()
    );
    println!("Inv Mass {} helix{}", vm.helices.len(), mm);

    let mm = inv_mass( &l5.iter()
                        .map( |i| PMeas::from(&QMeas::from(&vm.helices[*i])) )
                        .collect()
                    );
    println!("Inv Mass {} helix{}", l5.len(), mm);
    // for p in &pl5 { println!("{}", p) };

    println!("Fitting Vertex --------------------");
    let Prong { fit_vertex: vf,
                fit_momenta: qs,
                fit_chi2s: cs,
                n_prong: np,
                measurements: _ms
                } = fit(&vm);
    println!("Fitted vertex -> {}", vf);

    for i in 0..np { println!("q chi2 ->{:6.1} {}", cs[i], qs[i]); }

    let mut pl: Vec<PMeas> = Vec::new();
    for q in &qs { pl.push(PMeas::from(q)); }
    println!("Inv Mass {} fit{}", np, inv_mass(&pl));

    let mut pl5: Vec<PMeas> = Vec::new();
    for &i in &l5 { pl5.push(PMeas::from(&qs[i])); }
    println!("Inv Mass {} fit{}", l5.len(), inv_mass(&pl5));

    println!("Refitting Vertex-----------------");
    let h5s = l5.iter().map( |i| vm.helices[*i].clone() ).collect();
    let vmp = VHMeas{ helices: h5s, ..vm };
    let Prong {fit_vertex: fv,
        fit_momenta: fqs,
        fit_chi2s: fcs,
        n_prong: fnp,
        measurements: _} = fit(&vmp);
    println!("Refitted vertex -> {}", fv);
    for i in 0..fnp { println!("q chi2 ->{:6.1} {}", fcs[i], fqs[i]); }

    let mm = inv_mass( &fqs
            .iter()
            .map( |q| PMeas::from(q) )
            .collect()
        );
    println!("Inv Mass {} refit{}", fnp, mm);

    println!("Final vertex -> {}", fv);
    println!("end of doFitTest------------------------------------------");

    let res = String::from("all good?");
    assert!( 1 == 1, "test failed with '{}'", res);

}

