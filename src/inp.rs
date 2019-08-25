

use crate::types::*;
use crate::cov::*;

pub fn h_slurp(ds: String) -> Option<VHMeas> {
    let ws = ds.split_whitespace().collect::<Vec<&str>>();
    println!("h_slurp len = {:?}", ws.len());

    // sometimes there is PU information at the front -- skip for now
    let npu: Option<usize> = if ws[0] == "PU_zpositions:" { 
        Some(ws[1].parse().unwrap())
    } else { None };
    println!("h_slurp PU = {:?}", npu);

    let varr = match npu {
        None    => ws.iter().map(|f| f.parse::<f64>().unwrap()).collect(),
        Some(n) => ws[(n+2)..].iter().map(|f| f.parse::<f64>().unwrap()).collect(),
    };
    h_slurpp(varr)
}

fn h_slurpp(inp: Vec<f64>) -> Option<VHMeas> {
    let v0: Vec3   = inp[..3].to_vec().into();       // initial vertex pos
    let cv0: Cov3  = inp[3..12].to_vec().into();     // cov matrix
    let v    = XMeas(v0, cv0);
    println!("h_slurp v = {:?}", v);
    let w2pt = inp[12];                    // how to calc pt from w; 1 in case of CMS
    let nt   = inp[13] as usize;           // number of helices to follow
      // f     = case w2pt of
      //             1.0 -> nxtH'        -- CMS case
      //             otherwise -> nxtH   -- Aleph case
      // hl    = mapMaybe (\i -> f w2pt (slice (i*30+14) (i*30+44) inp)) $ range 0 (nt-1)

    println!("h_slurp w2pt = {:?}", w2pt);
    println!("h_slurp nt = {:?}", nt);
    let mut hl: Vec<HMeas> = Vec::new();
    let aleph = w2pt != 1.0;
    for i in 0..nt {
        let h0 =    if aleph { nxt_h(w2pt, inp[i*30+14..i*30+44].to_vec()).unwrap() }
                    else { nxt_hp(inp[i*30+14..i*30+44].to_vec()).unwrap() };
        hl.push(h0);
    }
    println!("h_slurp h0 = {:?}", hl[0]);
    Some(VHMeas{ vertex: v, helices: hl })
}

// -- get the next helix, aleph case
fn nxt_h(w0: Number, ds: Vec<Number>) -> Option<HMeas> {
    let h = ds[..5].to_vec().into();
    let ch = ds[5..30].to_vec().into();
    Some(HMeas(h, ch, w0))
}

// -- get the next helix, CMS case
fn nxt_hp(ds: Vec<Number>) -> Option<HMeas> {
  // -- FV works in terms of a perigee system
  // -- w = omega = 1/R is curvature radius
  // -- tl = tan lambda = tangent of dipping angle (0 for pt-max)
  // -- psi = angle in xy plane
  // -- d0 = distance of closest approach to 0 in xy plan
  // -- z0 = positiion on z axis
  // --
  // -- CMS uses instead
  // -- q/p = charge over momentum
  // -- theta = dip angle
  // -- etc
    let h0 = ds[0];
    let h1 = ds[1];
    let h2 = ds[2];
    let h3 = ds[3];
    let h4 = ds[4];
    let w0            = 0.003*3.8;  // CMS case: field is 3.8 T, give R in cm
    let st            = h1.sin();
    let ct            = h1.cos();
    let w             = h0 * w0 / ct;
    let tl            = st / ct;
    let j00           = w0 / ct;
    let j01           = h0 * w0 * st/ct/ct;
    let j11           = 1.0 / ct / ct;
    let j10           = 0.0;
    let jj: Jac55     = [ j00, j01, 0.0, 0.0, 0.0,
                          j10, j11, 0.0, 0.0, 0.0,
                          0.0, 0.0, 1.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 1.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 1.0,
                        ].into();
    let hp: Vec5        = [w, tl, h2, h3, h4].into();
    let chp: Cov5       = ds[5..30].into();
    let chpp            = jj * chp;

    Some(HMeas(hp, chpp, w0))

}


#[test]
fn test_inp_aleph() {
    let ds = std::fs::read_to_string("dat/tr05129e001412.dat").unwrap();
    let VHMeas {vertex: _x, helices: hl} = h_slurp(ds).unwrap();
    let HMeas(_x,_y, w) = &hl[hl.len()-1];

    let res = String::from("all good?");
    assert!( *w == 4.5451703e-3, "test failed with '{}'", res);
}
#[test]
fn test_inp_cms() {
    let _ds = std::fs::read_to_string("dat/tav-1.dat").unwrap();
    let ds = TAV4.to_string();
    let VHMeas {vertex: _x, helices: hl} = h_slurp(ds).unwrap();
    let HMeas(_x,_y, w) = &hl[hl.len()-1];

    let res = String::from("all good?");
    assert!( *w == 0.0114f64, "test failed with '{}'", res);
}

const TAV4: &'static str = r"PU_zpositions:  190 4.06972837448 2.44204807281 7.82136058807 -0.621172726154 -6.80061435699 -1.73116350174 -5.42739343643 -7.10662841797 -6.32562208176 -3.72315001488 1.66695046425 6.55822181702 -7.12538957596 -0.389555871487 -2.8334877491 3.09819436073 -5.65534687042 12.068236351 -1.79448211193 5.73383188248 1.68428444862 2.1804420948 8.66328144073 -12.8040647507 -1.1730145216 -3.57441878319 6.21948480606 -1.26211774349 -3.4871032238 -9.48501300812 -8.33902263641 -1.71619582176 -1.56027853489 1.49686825275 -1.69698286057 1.69038307667 5.10251283646 -2.57128977776 0.749759852886 -2.58463263512 -9.792719841 -8.84095287323 -0.131224393845 -1.56865620613 -5.81232976913 4.21827507019 -4.92665529251 -5.84215211868 -5.74135446548 3.38353490829 -3.13945651054 4.30185222626 -12.6121692657 1.54116880894 1.38944470882 -6.84423398972 2.88845825195 -4.16181087494 6.3093957901 -1.70226609707 3.62256598473 -1.38095474243 1.69552695751 -9.44017601013 2.82410240173 -2.21053552628 2.34878325462 -8.67048835754 1.25067412853 9.49777984619 8.16330623627 -0.870663702488 -4.79498910904 1.78941035271 -7.03154611588 1.68979644775 -0.484967201948 -4.18258905411 0.0788396298885 -4.69477128983 2.32463097572 -2.10498857498 -5.34199571609 3.32180857658 -5.39752531052 -2.84948658943 -2.68618583679 1.0778503418 0.443690419197 -3.29635429382 0.936188876629 -4.41851854324 -3.29131436348 2.12316703796 -10.6452322006 -14.0393047333 3.74121594429 -8.4497051239 -5.68886137009 8.31489753723 -4.49255418777 -7.92309999466 -7.26154613495 -2.43943715096 2.87128973007 -8.41958713531 -5.04697036743 -2.6269865036 -3.01578998566 5.666908741 4.7386713028 4.83959341049 -12.2599534988 6.80844593048 -7.59651374817 1.77152347565 -3.49425053596 4.14569759369 2.39712738991 0.695241510868 0.351206511259 -1.00542604923 -0.592145264149 8.05185890198 1.35937333107 -3.23685288429 1.82836604118 -1.08040130138 -4.06748771667 -1.22976350784 -5.24559354782 4.77764129639 -7.92655897141 6.87241268158 8.90295886993 -10.4462614059 5.51054620743 4.28739690781 -0.413518726826 -2.84266161919 -4.82323074341 -3.47484374046 -6.56179046631 -5.6174902916 2.68036007881 -4.87207984924 -3.47317409515 -1.94823920727 -11.0047950745 -6.04952716827 -12.1523780823 -0.171474739909 1.82068359852 -11.1572389603 -2.97859430313 -3.65392804146 1.67614769936 -4.62239599228 4.72258663177 -3.13622426987 -9.94389533997 -13.6851511002 1.98555517197 4.60026597977 -10.9611978531 -1.63044011593 8.50263690948 -9.76078033447 0.933302462101 6.68330335617 -2.94098043442 -8.59897899628 -0.908704698086 -5.6248884201 -9.19552707672 -6.67034435272 3.34288668633 -2.66896915436 -5.85388660431 -6.08788156509 -9.28157234192 -3.39719057083 -2.08446788788 3.61256814003 4.3055267334 -3.20882606506 -1.37032854557 6.3657708168 -7.99672412872 7.93814659119
0.104794 0.168646 -1.00377 0.0015033299569 0.0 0.0 0.0 0.00151841994375 0.0 0.0 0.0 5.21037006378
1.0
4
-0.450663641447 1.35035226203 1.18063795337 0.0660153061812 -1.23642665653 1.75648085587e-06 1.09397257919e-09 2.97465732046e-07 -4.61079963543e-07 -2.34128183507e-08 1.09397257919e-09 1.44003013247e-06 8.51150190329e-08 2.84054323174e-07 -1.69398772414e-05 2.97465732046e-07 8.51150190329e-08 3.07370719383e-05 -7.83819778007e-05 -3.34401283908e-06 -4.61079963543e-07 2.84054323174e-07 -7.83819778007e-05 0.000201181290322 2.55465602095e-06 -2.34128183507e-08 -1.69398772414e-05 -3.34401283908e-06 2.55465602095e-06 0.000202862953302
0.425837572723 -1.33656916571 -0.200853553729 0.186440212098 0.447653308602 4.94152754982e-06 -1.85234905192e-09 7.72127236814e-07 -1.24222322029e-06 -3.99287580777e-09 -1.85234905192e-09 1.2965380165e-06 1.72925282982e-08 4.64915046905e-07 -1.68021942955e-05 7.72127236814e-07 1.72925282982e-08 2.40928802668e-05 -7.24880374037e-05 -2.46094759859e-06 -1.24222322029e-06 4.64915046905e-07 -7.24880374037e-05 0.000219607783947 7.73318163283e-07 -3.99287580777e-09 -1.68021942955e-05 -2.46094759859e-06 7.73318163283e-07 0.000218317465624
0.292514034965 1.39665579944 -0.993975573833 0.141543926193 0.40182813437 6.18352771653e-07 2.45095588269e-09 2.34424803125e-07 -8.39045242174e-07 -8.35598825688e-08 2.45095588269e-09 8.66051891535e-07 4.61750582215e-10 -3.98558910319e-07 -1.4314922737e-05 2.34424803125e-07 4.61750582215e-10 2.8809732612e-05 -8.27241237857e-05 2.18257969209e-06 -8.39045242174e-07 -3.98558910319e-07 -8.27241237857e-05 0.000238582899328 4.61562564169e-07 -8.35598825688e-08 -1.4314922737e-05 2.18257969209e-06 4.61562564169e-07 0.000236972875427
-0.29652562498 1.37864870921 -1.02387889924 0.192178996941 0.0749212341852 1.4916007558e-06 5.99006577673e-09 7.09483913397e-07 -2.9230166092e-06 -3.92529898363e-07 5.99006577673e-09 7.67293840909e-07 5.73736658183e-09 2.58630308281e-07 -1.14107451736e-05 7.09483913397e-07 5.73736658183e-09 2.14707033592e-05 -6.15742173977e-05 -1.75877414677e-06 -2.9230166092e-06 2.58630308281e-07 -6.15742173977e-05 0.000178089467227 1.27894770685e-06 -3.92529898363e-07 -1.14107451736e-05 -1.75877414677e-06 1.27894770685e-06 0.000170099316165
";

