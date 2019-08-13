
use crate::*;
use crate::chol::*;

pub type NA3 = [Number; 3];
pub type NA6 = [Number; 6];
pub type NA9 = [Number; 9];

#[derive(Debug, PartialEq, Clone)]
pub struct Vec3 ( pub NA3 );

#[derive(Debug, PartialEq, Clone)]
pub struct Cov3 { pub v: NA6 }

impl Cov3 {
    pub fn to_string(&self) -> String {
        let w = 3;
        let idx = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let mut v: NA9 = [0.0; 9];
        for i in 0usize..3 {
            for j in 0usize..3 {
                v[j*w+i] = self.v[idx(i,j)];
            }
        }
        pretty_matrix(3,3,v)
    }

    // SymMat
    pub fn det(&self) -> Number {
        // [a,b,c,d,e,f] = self.v;
        let a = self.v[0];
        let b = self.v[1];
        let c = self.v[2];
        let d = self.v[3];
        let e = self.v[4];
        let f = self.v[5];
        a*d*f - a*e*e - b*b*f + 2.0*b*c*e - c*c*d
    }

    pub fn choldc(self) -> Jac33 {
        let vj: NA9 = do_choldc(self.v, 3);
        Jac33 { v: vj }
    }

    pub fn cholinv(&self) -> Cov3 {
        let vc: NA6 = do_cholinv(self.v, 3);
        Cov3 { v: vc }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Jac33 { pub v: NA9 }

impl Jac33 {
    pub fn to_string(&self) -> String {
        // format!("Jac33 {:?}", self.v)
        pretty_matrix(3, 3, self.v)
    }
    pub fn tr(mut self) -> Jac33 {

        let w = 3;
        let ixa = |i0: usize, j0: usize| i0*w+j0;

        let tmp: NA9 = self.v.clone();                ;
        for i0 in 0..w {
            for j0 in 0..w {
                self.v[ixa(i0, j0)] = tmp[ixa(j0, i0)];
            }
        }
        self
    }
}

// fn new_na9_from<F: Iterator<Item=Number>>(src: F) -> NA9 {
//     let mut result: NA9 = [0.0; 9];
//     for (rref, val) in result.iter_mut().zip(src) {
//         *rref = val;
//     }
//     result
// }
// fn cmap<F: Fn(Number, Number) -> Number>(a: NA9, b: NA9, f: F) -> NA9 {
//     new_na9_from(a.iter().zip(&b).map(|(x, y)| f(*x, *y)))
// }

use std::ops::Add;
impl Add for Cov3 {
    type Output = Cov3;
    fn add(mut self, other: Cov3) -> Cov3 {
        for i in 0..6 { self.v[i] += other.v[i]; }
        self
    }
}
impl Add for Jac33 {
    type Output = Jac33;
    fn add(mut self, other: Jac33) -> Jac33 {
        for i in 0..9 { self.v[i] += other.v[i]; }
        // Jac33 {
        //     // v: cmap(self.v, other.v, |x, y| x + y),
        //     // v: self.v.iter().zip(&other.v).map(|(s, o)| s+o).collect()
        // }
        self
    }
}
use std::ops::Mul;
impl Mul for Jac33 {
    type Output = Jac33;
    fn mul(mut self, other: Jac33) -> Jac33 {
        //self.v.len() is 9;
        let nb = 3;
        let na = self.v.len() / nb;
// indV w i0 j0 = i0*w+j0 -- w=nj width of niXnj matrix, i0=0..ni-1, j0=0..nj-1
        let ixa = |i0, j0| i0*nb+j0;
        let ixb = |i0, j0| i0*na+j0;
        let tmp = self.v.clone();                ;
        for i in 0..na {
        for j in 0..na {
        let mut s = 0.0;
        for k in 0..nb {
            s += tmp[ixa(i,k)] * other.v[ixb(k,j)];
        }
        self.v[ixa(i,j)] = s;
        }}
        self
    }
}

impl Mul for Cov3 {
    type Output = Jac33;
    fn mul(self, other: Cov3) -> Jac33 {
        //self.v.len() is 9;
        let n = 3;
        let w = n;
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let ixr = |i0: usize, j0: usize| i0*w+j0;

        let mut res = Jac33 { v: [0.0; 9] };
        for i in 0..n {
        for j in 0..n {
        let mut s = 0.0;
        for k in 0..n {
            s += self.v[ixa(i,k)] * other.v[ixa(k,j)];
        }
        res.v[ixr(i,j)] = s;
        }}
        res
    }
}

#[derive(Debug)]
pub struct Cov5(pub [Number; 15]);
impl Cov5 {
    pub fn to_string(&self) -> String {
        format!("Cov5 {:?}", self.0)
    }
    // SymMat
    pub fn det(&self) -> Number {
        // [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o] = self.0;
        let (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) = (self.0[0],self.0[1],self.0[2],self.0[3],self.0[4],self.0[5],self.0[6],self.0[7],self.0[8],self.0[9],self.0[10],self.0[11],self.0[12],self.0[13],self.0[14],);
        a*f*j*m*o - a*f*j*n*n - a*f*k*k*o + 2.0*a*f*k*l*n - a*f*l*l*m
            - a*g*g*m*o + a*g*g*n*n + 2.0*a*g*h*k*o - 2.0*a*g*h*l*n - 2.0*a*g*i*k*n
            + 2.0*a*g*i*l*m - a*h*h*j*o + a*h*h*l*l + 2.0*a*h*i*j*n - 2.0*a*h*i*k*l
            - a*i*i*j*m + a*i*i*k*k - b*b*j*m*o + b*b*j*n*n + b*b*k*k*o
            - 2.0*b*b*k*l*n + b*b*l*l*m + 2.0*b*c*g*m*o - 2.0*b*c*g*n*n - 2.0*b*c*h*k*o
            + 2.0*b*c*h*l*n + 2.0*b*c*i*k*n - 2.0*b*c*i*l*m - 2.0*b*d*g*k*o
            + 2.0*b*d*g*l*n + 2.0*b*d*h*j*o - 2.0*b*d*h*l*l - 2.0*b*d*i*j*n
            + 2.0*b*d*i*k*l + 2.0*b*e*g*k*n - 2.0*b*e*g*l*m - 2.0*b*e*h*j*n
            + 2.0*b*e*h*k*l + 2.0*b*e*i*j*m - 2.0*b*e*i*k*k - c*c*f*m*o + c*c*f*n*n
            + c*c*h*h*o - 2.0*c*c*h*i*n + c*c*i*i*m + 2.0*c*d*f*k*o - 2.0*c*d*f*l*n
            - 2.0*c*d*g*h*o + 2.0*c*d*g*i*n + 2.0*c*d*h*i*l - 2.0*c*d*i*i*k
            - 2.0*c*e*f*k*n + 2.0*c*e*f*l*m + 2.0*c*e*g*h*n - 2.0*c*e*g*i*m
            - 2.0*c*e*h*h*l + 2.0*c*e*h*i*k - d*d*f*j*o + d*d*f*l*l + d*d*g*g*o
            - 2.0*d*d*g*i*l + d*d*i*i*j + 2.0*d*e*f*j*n - 2.0*d*e*f*k*l - 2.0*d*e*g*g*n
            + 2.0*d*e*g*h*l + 2.0*d*e*g*i*k - 2.0*d*e*h*i*j - e*e*f*j*m + e*e*f*k*k
            + e*e*g*g*m - 2.0*e*e*g*h*k + e*e*h*h*j
    }
}


/// pretty print of matrix
pub fn pretty_matrix(r: usize, c: usize, v: NA9) -> String {
    let to3fix = |x: Number| format!("{:.3}", x);
    let fill_blanks = |k: usize, str: &str| format!("{:>1$}", str, k);
    let mx: usize = v.iter().map(|x| to3fix(*x).len()).max().unwrap();
    let idx = |i0, j0| j0*c+i0;
    let fmt = |x: Number| fill_blanks(mx, &to3fix(x));
    (0usize..r).map(|j| format!("({})",
                                (0usize..c).
                                map(|i| fmt(v[idx(i,j)])).
                                fold("".to_string(), |accum, s| accum + &s + " ")
                                )
                    ).fold("".to_string(), |a, s| a + "\n" + &s)
}

#[test]
fn test_cov() {
    // let Cov<Dim3> xc3 = Cov {v: Vec}
    //
    let ch3 = Cov3 {v: [2.0, -1.0, 0.0, 2.0, -1.0, 2.0]};
    let ch5 = Cov5([2.0, -1.0, 0.0, 0.0, 0.0, 2.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0]);


    let res = format!(r"
chol: -----------------
A = L * L^T             {}
L                       {}
L * L^T                 {}
A^(-1) = L' * L'^T      {}
A * A^(-1)              {}
A = L * L^T             {}
L                       {}
L * L^T                 {}
A^(-1) = L' * L'^T      {}
A * A^(-1)              {}
det this                {}
",
    ch3.to_string(),
    // ch3.to_string(),
    ch3.clone().choldc().to_string(),
    // (choldc ch3).to_string(),
    (ch3.clone().choldc() * ch3.clone().choldc().tr()).to_string(),
    // ((choldc ch3) *. tr (choldc ch3)).to_string(),
    ch3.cholinv().to_string(),
    // (cholInv ch3).to_string(),
    (ch3.clone() * ch3.cholinv()).to_string(),
    // (ch3 *. cholInv ch3).to_string(),
    ch5.to_string(),
    // ch5.to_string(),
    ch3.to_string(),
    // (choldc ch5).to_string(),
    ch3.to_string(),
    // ((choldc ch5) *. tr (choldc ch5)).to_string(),
    ch3.to_string(),
    // (cholInv ch5).to_string(),
    ch3.to_string(),
    // (ch5 *. cholInv ch5).to_string(),
    ch5.det().to_string(),
    // ch5.det().to_string(),
    );
    print!("{}", res);
    assert!(true, "test failed with '{}'", res);
}


