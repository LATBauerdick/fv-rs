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

use crate::inp::h_slurp;
#[test]
fn test_fvt() {
    println!("test_fvt-------------------------------------------------");
    let ds = std::fs::read_to_string("dat/tr05129e001412.dat").unwrap();
    let VHMeas {vertex: x, helices: hel} = h_slurp(ds).unwrap();
    for h in &hel { println!("{}", h) };
    for h in &hel { println!("{}", QMeas::from(h)) };
//   doFitTest vm l5
//   putStrLn $ showProng <<< fit <<< hFilter l5 <<< vBlowup 10000.0 $ vm


    println!("init vtx pos -> {}", x);

    let vm = VHMeas {vertex: x.blowup(10000.0), helices: hel};
    // let pl = map(fromQMeas <<< fromHMeas) $ helices


    println!("---------------------------------------------------------");
    let res = String::from("all good?");
    assert!( 1 == 1, "test failed with '{}'", res);

}

// doFitTest :: VHMeas
//             -> List Int
//             -> IO ()
// doFitTest vm' l5 = do
//   let vm = vBlowup 10000.0 vm'
//   let showLen xs = show $ length xs
//       showQChi2 :: (QMeas, Chi2) -> String
//       showQChi2 (qm, (Chi2 chi2)) = "q"
//                                 <> " chi2 ->" <> to1fix chi2
//                                 <> " pt,pz,fi,E ->"
//                                 <> show qm

//   putStrLn $           "initial vertex position -> " <> show ((vertex vm)::XMeas)

//   let pl         = map (fromQMeas <<< fromHMeas) $ helices vm
//   putStrLn $ "Inv Mass " <> showLen pl <> " helix" <> show (invMass pl)
//   let pl5        = map (fromQMeas <<< fromHMeas) (helices <<< hFilter l5 $ vm)
//   putStrLn $ "Inv Mass " <> showLen pl5 <> " helix" <> show (invMass pl5)

//   putStrLn             "Fitting Vertex --------------------"
//   let -- pr = fit vm
//       Prong {fitVertex= vf, fitMomenta= ql, fitChi2s= cl} = fit vm
//   putStrLn $           "Fitted vertex -> " <> show vf
//   traverse_ (putStrLn <<< showQChi2) $ zip ql cl
//   putStrLn $ "Inv Mass " <> show (length ql) <> " fit"
//                     <> show (invMass (map fromQMeas ql))

//   let m5 = invMass <<< map fromQMeas <<< iflt l5 $ ql
//   putStrLn $ "Inv Mass " <> show (length l5) <> " fit" <> show m5

//   putStrLn $           "Refitting Vertex-----------------"
//   let Prong {fitVertex=fv, fitMomenta=fqs, fitChi2s=fcs, nProng=np} = fit <<< hFilter l5 $ vm
//   putStrLn $           "Refitted vertex -> " <> show fv
//   traverse_ (putStrLn <<< showQChi2) $ zip fqs fcs
//   putStrLn $           "Inv Mass " <> show np <> " refit" 
//                        <> (show <<< invMass <<< map fromQMeas $ fqs)
//   putStrLn $           "Final vertex -> " <> show fv
//   putStrLn $           "end of doFitTest------------------------------------------"
