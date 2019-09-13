
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

#[derive(Debug, Clone)]
pub struct VHMeas {
    pub vertex:  XMeas,
    pub helices: Vec<HMeas>,
}

#[derive(Debug, Clone)]
pub struct Chi2(pub Number);
impl fmt::Display for Chi2 {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        write!(fmt, "{:6.1}", self.0)
    }
}

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
// -- | calculate q 3-vector for a given helix parameterization near vertex position
    pub fn hv2q(h: &Vec5, v: &Vec3) -> Vec3 {
        let xx   = v.v[0];
        let yy   = v.v[1];
        let r    = f64::sqrt(xx*xx + yy*yy);
        let phi  = f64::atan2(yy, xx);
        let w0   = h.v[0];
        let tl0  = h.v[1];
        let psi0 = h.v[2];
//   -- let d0   = h.v[2];
//   -- let z0   = h.v[2];
        let xi   = ((psi0 - phi) % TWOPI + TWOPI) % TWOPI;
        let (sxi, cxi)   = f64::sin_cos(xi);
        let qv = if w0 != 0.0 {
                                let oow0 = 1.0/w0;
                                let gamma = f64::atan(r*cxi/(oow0-r*sxi));
                                [ w0, tl0, psi0 + gamma ]
                            } else {
                                [ w0, tl0, psi0 ]
                            };

        // let qv = [h.v[0],h.v[1],h.v[2]];
        Vec3 { v: qv }
    }
// hv2q Data.Cov.Vec.Vec {Data.Cov.Vec.v=h_} Data.Cov.Vec.Vec {Data.Cov.Vec.v=v_} = q where


// hv2q :: Vec5 -> Vec3 -> Vec3
// hv2q Data.Cov.Vec.Vec {Data.Cov.Vec.v=h_} Data.Cov.Vec.Vec {Data.Cov.Vec.v=v_} = q where


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
        let cqp        = &jj % &cq;
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
    fn mass(&self) -> MMeas {
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
        PMeas(&self.0+&other.0, &self.1+&other.1)
    }
}
pub fn inv_mass(ps: &Vec<PMeas>) -> MMeas {
    let psum = ps[1..].iter()
                      .fold( ps[0].clone(), |acc, p| acc + p ); // PMeas should get a Default
    psum.mass()
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

const TWOPI: f64 = 2.0*PI;
pub fn expand(v: &Vec3, q: &Vec3) -> ( Jac53, Jac53, Vec5 ) {
    let xx  = v.v[0];
    let yy  = v.v[1];
    let z   = v.v[2];
    let r   = f64::sqrt(xx*xx + yy*yy);
    let phi = f64::atan2(yy, xx);
    let w   = q.v[0];
    let tl  = q.v[1];
    let psi = q.v[2];
  // -- some more derived quantities
    let xi  = ((psi - phi) % TWOPI + TWOPI) % TWOPI;
    let (sxi, cxi)   = f64::sin_cos(xi);
    let oow = 1.0 / w;
    let rw  = r * w;

    let gamma = f64::atan(r*cxi/(oow - r*sxi));
    let oow = 1.0 / w;
    let (sg, cg)   = f64::sin_cos(gamma);

  // -- calculate transformed quantities
    let psi0  = psi - gamma;
    let d0    = oow - (oow - r*sxi)/cg;
    let z0    = z - tl*gamma/w;

  // -- calc Jacobian
    let drdx    =    if r != 0.0 {  xx/r } else { 0.0 };
    let drdy    =    if r != 0.0 {  yy/r } else { 0.0 };
    let rdxidx  =    if r != 0.0 {  yy/r } else { 0.0 };
    let rdxidy  =    if r != 0.0 { -xx/r } else { 0.0 };
    let dgdvar0 =    1.0/(1.0 + rw*rw - 2.0*rw*sxi);
    let dgdx    =    dgdvar0*(w*cxi*drdx + w*(rw - sxi)*rdxidx);
    let dgdy    =    dgdvar0*(w*cxi*drdy + w*(rw - sxi)*rdxidy);
    let dgdw    =    dgdvar0*r*cxi;
    let dgdpsi  =    dgdvar0*rw*(rw - sxi);

  // --  fill matrix:
  // -- d w / d r, d phi, d z
    let a11                = 0.0;
    let a12                = 0.0;
    let a13                = 0.0;
  // -- d tl / d x, d y, d z
    let a21                = 0.0;
    let a22                = 0.0;
    let a23                = 0.0;
  // -- d psi0 / d x, d y, d z
    let a31                = -dgdx;
    let a32                = -dgdy;
    let a33                = 0.0;
  // -- d d0 / d x, d y, d z
    let a41                = cxi*rdxidx/cg + sxi*drdx/cg
                                - (oow - r*sxi)*sg*dgdx/cg/cg;
    let a42                = cxi*rdxidy/cg + sxi*drdy/cg
                                - (oow - r*sxi)*sg*dgdy/cg/cg;
    let a43                = 0.0;
  // -- d z0 / d x, d y, d z
    let a51                = -tl/w*dgdx;
    let a52                = -tl/w*dgdy;
    let a53                = 1.0;

  // -- B
  // -- d w / d w, d tl, d psi
    let b11                = 1.0;
    let b12                = 0.0;
    let b13                = 0.0;
  // -- d tl / d w, d tl, d psi
    let b21                = 0.0;
    let b22                = 1.0;
    let b23                = 0.0;
  // -- d psi0 / d w, d tl, d psi
    let b31                = -dgdw;
    let b32                = 0.0;
    let b33                = 1.0 - dgdpsi;
  // -- d d0 / d w, d tl, d psi
    let b41                =  -oow*oow*(1.0 - 1.0/cg)
                       - (oow - r*sxi)*sg*dgdw/cg/cg;
    let b42                = 0.0;
    let b43                = r*cxi/cg - (oow - r*sxi)*sg*dgdpsi/cg/cg;
  // -- d z0 / d w, d tl, d psi
    let b51                = -tl/w*(dgdw - gamma/w);
    let b52                = -gamma/w;
    let b53                = -tl/w*dgdpsi;

    let  v01                = xx;
    let  v02                = yy;
    let  v03                = z;
    let  q01                = w;
    let  q02                = tl;
    let  q03                = psi;
    let  h0                 = Vec5 { v: [
                        0.0,
                        0.0,
                        psi0 - a31*v01 - a32*v02 - b31*q01 - b33*q03,
                        d0 - a41*v01 - a42*v02 - b41*q01 - b43*q03,
                        z0 - a51*v01 - a52*v02 - a53*v03 - b51*q01 - b52*q02 - b53*q03
                    ] };
    let  aa = Jac53 { v: [a11,a12,a13,a21,a22,a23,a31,a32,a33,a41,a42,a43,a51,a52,a53] };
    let  bb = Jac53 { v: [b11,b12,b13,b21,b22,b23,b31,b32,b33,b41,b42,b43,b51,b52,b53] };

    // let aa= Jac53::default();
    // let bb= Jac53::default();
    // let h0 = Vec5::default();
    (aa, bb, h0)
}







