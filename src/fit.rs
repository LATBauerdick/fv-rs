

use crate::cov::*;
use crate::types::*;

use std::fmt;

use std::convert::From;

pub fn fit<'a>(vhm: &'a VHMeas) -> Prong<'a> {
    vhm.k_smooth(vhm.k_filter())
}

impl VHMeas {
    fn k_filter(&self) -> XMeas { self.vertex.clone() }
// k_filter VHMeas {vertex=v, helices=hl} = foldl k_add v hl

    // fn k_add() -> XMeas {

    // }
// kAdd :: XMeas -> HMeas -> XMeas
// kAdd (XMeas v vv) (HMeas h hh w0) = kAdd' x_km1 p_k x_e q_e (Chi2 1e6) 0 where
//   x_km1 = XMeas v (inv vv)
//   p_k   = HMeas h (inv hh) w0
//   x_e   = v
//   q_e   = J.hv2q h v

    fn k_smooth(&self, v: XMeas) -> Prong {
        let n = self.helices.len();
        let mut ql: Vec<QMeas> = Vec::new();
        for i in 0..n { ql.push(QMeas::from(&self.helices[i])); }
        Prong {n_prong: 0,
        fit_vertex: v,
        fit_momenta: ql,
        fit_chi2s: vec![1.0; n],
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

