
use crate::cov::*;

use std::fmt;

use std::convert::From;

pub type Number = f64;

#[derive(Debug, Clone)]
pub struct Prong<'a> {
                    pub n_prong: usize,
                    pub fit_vertex: XMeas,
                    pub fit_momenta: Vec<QMeas>,
                    pub fit_chi2s: Vec<Chi2>,
                    pub measurements: &'a VHMeas,
                }

type Chi2 = Number;

#[derive(Debug, Clone)]
pub struct XMeas(pub Vec3, pub Cov3);
impl XMeas {
    pub fn blowup(&self, scale: f64) -> XMeas {
        XMeas(self.0.clone(), self.1.scale_diag(scale))
    }
}
impl fmt::Display for XMeas {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let XMeas(v, cv) = self;
// -- return a string showing vertex position vector with errors

        let vv         = v.v;
        let x          = vv[0];
        let y          = vv[1];
        let z          = vv[2];
        let s2v: Vec<Number> = cv.diag().to_vec().into_iter().map(|x| x.sqrt()).collect();
        let dx         = s2v[0];
        let dy         = s2v[1];
        let dz         = s2v[2];
        fn f(x: &Number, dx: &Number) -> String {
            format!("{:7.2} +-{:7.2}", *x, *dx)
        }
        write!(fmt, "(r,z) =({:7.2}, {:7.2}), x y z ={}{}{}", f64::sqrt(x*x + y*y), z, f(&x, &dx), f(&y, &dy), f(&z, &dz))
    }
}

#[derive(Debug, Clone)]
pub struct HMeas(pub Vec5, pub Cov5, pub Number);
impl HMeas {

}
impl fmt::Display for HMeas {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let HMeas(h, ch, _w) = self;

        let hs = h.v;
        let sh: Vec<Number> = ch.diag().to_vec().into_iter().map(|x| x.sqrt()).collect();
        let s00 = format!("{:10.5} +-{:10.5}", hs[0], sh[0]);
        let s01 = format!("{:8.3} +-{:8.3}", hs[1], sh[1]);
        let s02 = format!("{:8.3} +-{:8.3}", hs[2], sh[2]);
        let s03 = format!("{:8.3} +-{:8.3}", hs[3], sh[3]);
        let s04 = format!("{:8.3} +-{:8.3}", hs[4], sh[4]);

        write!(fmt, "Helix ->{}{}{}{}{}", s00, s01, s02, s03, s04)
    }
}

#[derive(Debug, Clone)]
pub struct QMeas(pub Vec3, pub Cov3, pub Number);
impl QMeas {

}
impl From<&HMeas> for QMeas {
    fn from(hm: &HMeas) -> Self {
        let h = &hm.0;
        let ch = &hm.1;
        let w = &hm.2;
        let cq:Cov3 = ch.into();
        let q: Vec3 = [h.v[0],h.v[1],h.v[2]].into();
        QMeas(q,cq,*w)
    }
}
static MPI: f64 = 0.1395675_f64;
use std::f64::consts::PI;
impl fmt::Display for QMeas {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let QMeas(q, cq, w2pt) = self;
        fn f(s: &String, (x, dx): (&Number, &Number)) -> String {
            format!("{}{:8.3} +-{:8.3}", s, *x, *dx)
        }
        let m           = MPI;
        let wp          = w2pt;
        let w           = q.v[0];
        let tl          = q.v[1];
        let psi0        = q.v[2];
        let pt          = wp / w.abs();
        let pz          = pt*tl;
        let psi         = psi0*180.0/PI;
        let e           = f64::sqrt(pt*pt  + pz*pz + m*m);
        let jj   = Jac34 { v : [ -wp/w/w, -wp/w/w*tl, 0.0, -(pz*pz + pt*pt)/w/e
                                , 0.0, wp/w, 0.0, pt*pt*tl/e
                                , 0.0, 0.0, 1.0, 0.0] };
        let cqp        = jj * cq.clone();
        let pp         = [pt, pz, psi, e];
        let dp: Vec<Number> = cqp.diag().to_vec().into_iter().map(|x| x.sqrt()).collect();
        let dpp        = [dp[0], dp[1], dp[2]*180.0/PI, dp[3]];
        let sp         = pp[..].iter().zip(&dpp).fold("".to_string(), |s, x|{ f(&s, x) });
        write!(fmt, "qt,qz,fi,E ->{} GeV", sp)
    }
}

#[derive(Debug, Clone)]
pub struct MMeas{pub m: Number, pub dm: Number,}
impl fmt::Display for MMeas {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let MMeas{m, dm} = self;
        write!(fmt, " {:6.1} +-{:6.1} MeV", m*1000.0, dm*1000.0)
    }
}

#[derive(Debug, Clone)]
pub struct PMeas(pub Vec4, pub Cov4);
impl PMeas {
    fn pmass(&self) -> MMeas {
    let p     = &self.0;
    let cp    = &self.1;
    let px    = p.v[0];
    let py    = p.v[1];
    let pz    = p.v[2];
    let e     = p.v[3];
    let c11   = cp.v[0];
    let c12   = cp.v[1];
    let c13   = cp.v[2];
    let c14   = cp.v[3];
    let c22   = cp.v[4];
    let c23   = cp.v[5];
    let c24   = cp.v[6];
    let c33   = cp.v[7];
    let c34   = cp.v[8];
    let c44   = cp.v[9];
    let m     = (e*e-px*px-py*py-pz*pz).max(0_f64).sqrt();
    let sigm0 = px*c11*px + py*c22*py + pz*c33*pz + e*c44*e +
            2.0*(px*(c12*py + c13*pz - c14*e)
               + py*(c23*pz - c24*e)
               - pz*c34*e);
    let dm    = sigm0.max(0_f64).sqrt() / m;
        MMeas{m, dm}
}

}
impl fmt::Display for PMeas {
// -- print PMeas as a 4-momentum vector px,py,pz,E with errors
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let PMeas(p, cp) = self;
        let sp: Vec<Number> = cp.diag().to_vec().into_iter().map(|x| x.sqrt()).collect();
        write!(fmt, "px,py,pz,E -> {:8.3} +-{:8.3}{:8.3} +-{:8.3}{:8.3} +-{:8.3}{:8.3} +-{:8.3}", p.v[0], sp[0], p.v[1], sp[1], p.v[2], sp[2], p.v[3], sp[3])
    }
}
use std::ops::Add;
impl Add<&PMeas> for PMeas {
    type Output = PMeas;
    fn add(self, other: &PMeas) -> PMeas {
        PMeas(self.0+&other.0, self.1+&other.1)
    }
}
pub fn inv_mass(pl: Vec<PMeas>) -> MMeas {
    let mut ps = pl[0].clone();
    for p in &pl[1..] { ps = ps + p; }
    ps.pmass()
}
use std::f64;
impl From<&QMeas> for PMeas {
    fn from(qm: &QMeas) -> Self {
        let q    = &qm.0;
        let cq   = &qm.1;
        let w2pt = &qm.2;
        let m     = MPI;
        let w     = q.v[0];
        let tl    = q.v[1];
        let psi0  = q.v[2];
        let (sph, cph)   = f64::sin_cos(psi0);
        let  pt   = w2pt / w.abs();
        let px   = pt * cph;
        let py   = pt * sph;
        let pz   = pt * tl;
        let sqr  = |x| x*x;
        let e    = f64::sqrt(sqr(px) + py*py + pz*pz + m*m);
        let ps   = w2pt / w;
        let dpdk = ps*ps/w2pt;
        let c11  = cq.v[0];
        let c12  = cq.v[1];
        let c13  = cq.v[2];
        let c22  = cq.v[3];
        let c23  = cq.v[4];
        let c33  = cq.v[5];
        let xy   = 2.0*ps*dpdk*cph*sph*c13;
        let sxx  = sqr (dpdk*cph) * c11 + sqr (ps*sph) * c33 + xy;
        let  sxy  = cph*sph*(dpdk*dpdk*c11 - ps*ps*c33) +
           ps*dpdk*(sph*sph-cph*cph)*c13;
        let  syy  = sqr (dpdk*sph) * c11 + sqr (ps*cph) * c33 - xy;
        let  sxz  = dpdk*dpdk*cph*tl*c11 -
           ps*dpdk*(cph*c12-sph*tl*c13) -
           ps*ps*sph*c23;
        let  syz  = dpdk*dpdk*sph*tl*c11 -
           ps*dpdk*(sph*c12 + cph*tl*c13) +
           ps*ps*cph*c23;
        let  szz  = sqr (dpdk*tl) * c11 + ps*ps*c22 -
           2.0*ps*dpdk*tl*c12;
        let  sxe  = (px*sxx + py*sxy + pz*sxz)/e;
        let  sye  = (px*sxy + py*syy + pz*syz)/e;
        let  sze  = (px*sxz + py*syz + pz*szz)/e;
        let  see  = (px*px*sxx + py*py*syy + pz*pz*szz +
         2.0*(px*(py*sxy + pz*sxz) + py*pz*syz))/e/e;
        let cp: Cov4  = [ sxx, sxy, sxz, sxe
                             , syy, syz, sye
                                  , szz, sze
                                       , see].into();
        let p: Vec4   = [px,py,pz,e].into();
        PMeas(p,cp)
    }
}
#[derive(Debug, Clone)]
pub struct VHMeas {
    pub vertex:  XMeas,
    pub helices: Vec<HMeas>,
}

