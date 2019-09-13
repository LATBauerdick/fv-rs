

use crate::cov::*;
use crate::types::*;

// use std::fmt;

// use std::convert::From;

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
// -- | add a helix measurement to kalman filter, return updated vertex position
// -- | if we can't invert, don't update vertex
    fn k_add( XMeas(v0, vv0): XMeas, HMeas(h, hh, _w0): &HMeas ) -> XMeas {
        let uu0        = &vv0.cholinv();
        let gg         = &hh.cholinv();
        let mut q_e    = HMeas::hv2q(h, &v0);
        let mut x_e    = v0.clone();
        let mut chi2_0 = 1e6_f64;
        let mut iter   = 0;
        loop {
            let (aa, bb, h0) = expand(&x_e, &q_e);
            let ww   = (&bb % gg).cholinv();
            let gb   = gg - &(gg % &(&bb % &ww));
            let uu   = uu0 + &(&aa % &gb);
            let cc   = uu.cholinv();
            let p    = h - &h0;
            let v    = &cc * &(&(uu0 * &v0) + &(&aa * &(&gb * &p)));
            let dp   = &p - &(&aa * &v);
            let q    = &ww * &(&bb * &(gg * &dp));
            let dh   = &dp - &(&bb * &q);
            let dv   = &v - &v0;
            let chi2 = &dh * &(gg * &dh) + &dv * &(uu0 * &dv);

            const CHI2CUT: f64 = 0.5;
            const ITERMAX: usize = 99;
            let good_enough = f64::abs(chi2 - chi2_0) < CHI2CUT || iter > ITERMAX;

            if good_enough { return XMeas(v, cc); }
            chi2_0 = chi2;
            iter = iter+1;
            x_e = v;
            q_e = q;
        }
    }

    fn k_smooth(&self, v: XMeas) -> Prong {
        let n = self.helices.len();
        let mut ql: Vec<QMeas> = Vec::new();
        let mut cl: Vec<Chi2>  = Vec::new();
        for i in 0..n { if let Some((q,c)) = VHMeas::ksm(&v, &self.helices[i]) { ql.push(q); cl.push(c); }}
        // for i in 0..n { ql.push(QMeas::from(&self.helices[i])); }
        Prong {n_prong: n,
        fit_vertex: v,
        fit_momenta: ql,
        fit_chi2s: cl,
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


    // -- kalman smoother step: calculate 3-mom q and chi2 at kalman filter'ed vertex
    // -- if we can't invert, return Nothing and this track will not be included
    fn ksm(XMeas(x, cc): &XMeas, HMeas(h, hh, w0): &HMeas) -> Option<(QMeas, Chi2)> {
        let q_e    = HMeas::hv2q(h, x);
        let (aa, bb, h0) = expand(x, &q_e);
        let gg         = &hh.cholinv();
        let ww         = (&bb % gg).cholinv();
        let p          = h - &h0;
        let uu         = cc.cholinv();
        let dp   = &p - &(&aa * x);
        let q    = &ww * &(&bb * &(gg * &dp));
        let e0: Jac35  = cc * &aa;
        let e1: Jac35  = &e0 * gg;
        let e2: Jac53  = &bb * &ww;
        let ee: Jac33  = &e1 * &e2;
        let dd   = &ww + &(&ee % &uu);
        let r    = &p - &(&(&aa * x) + &(&bb * &q));
        let ch   = &r * &(gg * &r);

        let gb   = gg - &(gg % &(&bb % &ww));
        let uup  = &uu - &(&aa % &gb);
        let ccp  = uup.cholinv();
        let xp   = &ccp * &( &(&uu * x) - &(&aa *&(&gb * &p)));
        let dx   = x - &xp;
        let cx   = &dx * &(&uup * &dx);
        let chi2 = &cx + &ch;

        Some((QMeas(q, dd, *w0), Chi2(chi2)))
    }

}
