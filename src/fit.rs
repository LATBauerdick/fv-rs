

use crate::cov::*;
use crate::types::*;

// use std::fmt;

use std::convert::From;

pub fn fit<'a>(vhm: &'a VHMeas) -> Prong<'a> {
    vhm.k_smooth(vhm.k_filter())
}

impl VHMeas {
    // fn k_filter(&self) -> XMeas { self.vertex.clone() }
    fn k_filter(&self) -> XMeas {
        self.helices
            .iter()
            .fold(self.vertex.clone(), |v, h| VHMeas::k_add(v, &h) )
    }
    fn k_add( XMeas(v, vv): XMeas, HMeas(h, hh, w0): &HMeas ) -> XMeas {
        let q_e   = HMeas::hv2q(h, &v);
        let x_e   = v.clone();
        let x_km1 = XMeas(v, vv.cholinv());
        let p_k   = &HMeas(h.clone(), hh.cholinv(), *w0);
        VHMeas::k_addp(x_km1, p_k, x_e, q_e, Chi2(1e6_f64), 0)
    }
// -- | add a helix measurement to kalman filter, return updated vertex position
// -- | if we can't invert, don't update vertex
    fn k_addp(XMeas(v0, uu0): XMeas,
              HMeas(h, gg, _w0): &HMeas,
              x_e0: Vec3,
              q_e0: Vec3,
              Chi2(chi2_00): Chi2,
              iter0: usize
              ) -> XMeas {
        let mut x_e = x_e0;
        let mut q_e = q_e0;
        let mut chi2_0 = chi2_00;
        let mut iter = iter0;
        loop
        {
            let (aa, bb, h0) = expand(&x_e, &q_e);
            let ww: Cov3   = (&bb * gg).cholinv();
            let gb: Cov5   = gg - &(gg * &(&bb * &ww));
            let uu   = &uu0 + &(&aa * &gb);
            let cc   = uu.cholinv();
            let m    = h - &h0;
            let v    = &cc * &(&(&uu0 * &v0) + &(&aa * &(&gb * &m)));
            let dm   = &m - &(&aa * &v);
            let q    = &ww * &(&bb * &(gg * &dm));
            let dh   = &dm - &(&bb * &q);
            let dv   = &v - &v0;
            let chi2 = &dh * &(gg * &dh) + &dv * &(&uu0 * &dv);

            const CHI2CUT: f64 = 0.5;
            const ITERMAX: usize = 99;
            let good_enough = f64::abs(chi2 - chi2_0) < CHI2CUT || iter > ITERMAX;

            if good_enough {
                return XMeas(v, cc);
            } else {
                chi2_0 = chi2;
                iter = iter+1; // +1;
                x_e = v;
                q_e = q;
                // k_addp(XMeas(v0, uu0), HMeas(*h, *gg, *w0), v, q, Chi2(chi2), iter+1)
            }
        }
        // XMeas(v0, uu0)
    }
// -- | add a helix measurement to kalman filter, return updated vertex position
// -- | if we can't invert, don't update vertex
// kAdd' :: XMeas -> HMeas -> Vec3 -> Vec3 -> Chi2 -> Int -> XMeas
// --kAdd' (XMeas v0 uu0) (HMeas h gg w0) x_e q_e _ i |
// --        i == 0 && trace ("kadd'-->" <> show i <> "|" <> show v0 <> show h) false = undefined
// kAdd' (XMeas v0 uu0) (HMeas h gg w0) x_e q_e (Chi2 𝜒2_0) iter
//   | goodEnough = XMeas v cc
//   | otherwise  = kAdd' (XMeas v0 uu0) (HMeas h gg w0) v q (Chi2 𝜒2) (iter+1)
//   where
//     Jacs {aajacs=aa, bbjacs=bb, h0jacs=h0} = J.expand x_e q_e
//     aaT  = tr aa
//     bbT  = tr bb
//     ww   = fromJust $ invMaybe (bb .*. gg)
//     gb   = gg - gg .*. (bbT .*. ww)
//     uu   = uu0 + aa .*. gb
//     cc   = inv uu
//     m    = h - h0
//     v    = cc *. (uu0 *. v0 + aaT *. gb *. m)
//     dm   = m - aa *. v
//     q    = ww *. (bbT *. gg *. dm)
//     𝜒2   = (dm - bb *. q) .*. gg + (v - v0) .*. uu0
//     goodEnough = abs (𝜒2 - 𝜒2_0) < chi2cut || iter > iterMax where
//       chi2cut = 0.5
//       iterMax = 99 :: Int


    fn k_smooth(&self, v: XMeas) -> Prong {
        let n = self.helices.len();
        let mut ql: Vec<QMeas> = Vec::new();
        for i in 0..n { ql.push(QMeas::from(&self.helices[i])); }
        Prong {n_prong: n,
        fit_vertex: v,
        fit_momenta: ql,
        fit_chi2s: vec![Chi2(1.0); n],
        measurements: &self,
        }
    }

// --kSmooth vm v | trace ("kSmooth " <> (show <<< length <<< helices $ vm) <> ", vertex at " <> (show v) ) false = undefined
// kSmooth (VHMeas {vertex= v0, helices= hl}) v = pr' where
//   (ql, chi2l) = unzip $ mapMaybe (ksm v) hl
//   hl' = hl
//   n = length hl
//   n' = length ql
//   n'' = if n == n' then n else n' `debug` "kSmooth killed helices"
//   pr' = Prong { fitVertex= v, fitMomenta= ql, fitChi2s= chi2l, nProng= n'', measurements= VHMeas {vertex= v0, helices= hl'} }
    // }

}

