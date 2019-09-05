

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
        let q_e   = &HMeas::hv2q(h, &v);
        let x_e   = &v.clone();
        let x_km1 = XMeas(v, vv.cholinv());
        let p_k   = &HMeas(h.clone(), hh.cholinv(), *w0);
        VHMeas::k_addp(x_km1, p_k, x_e, q_e, Chi2(1e6_f64), 0)
    }
// -- | add a helix measurement to kalman filter, return updated vertex position
// -- | if we can't invert, don't update vertex
    fn k_addp(v0: XMeas, h: &HMeas, x_e: &Vec3, q_e: &Vec3, _: Chi2, i: usize) -> XMeas {
        v0
    }


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

