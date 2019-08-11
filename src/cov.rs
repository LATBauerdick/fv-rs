
use crate::*;
use crate::chol::*;

#[derive(Debug)]
pub struct Vec3 ( pub [Number; 3] );

#[derive(Debug)]
pub struct Cov3 {
    pub v: [Number; 6],
}

impl Cov3 {
    pub fn to_string(&self) -> String {
        format!("Cov3 {:?}", self.v)
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

    pub fn choldc(&self) -> Jac33 {
        let vj: &mut [Number; 9] = &mut [0.0; 9];
        do_choldc(self.v, 3, vj);
        Jac33 { v: *vj }
    }
}

#[derive(Debug)]
pub struct Jac33 {
    pub v: [Number; 9],
}
impl Jac33 {
    pub fn to_string(&self) -> String {
        format!("Jac33 {:?}", self.v)
    }
    pub fn tr(&mut self) -> &Jac33 {

        let w = 3;
        let ixa = |i0: usize, j0: usize| i0*w+j0;

        let tmp: [Number; 9] = self.v.clone();                ;
        for i0 in 0..w {
            for j0 in 0..w {
                self.v[ixa(i0, j0)] = tmp[ixa(j0, i0)];
            }
        }
        self
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
    ch3.choldc().to_string(),
    // (choldc ch3).to_string(),
    ch3.choldc().tr().to_string(),
    // ((choldc ch3) *. tr (choldc ch3)).to_string(),
    ch3.to_string(),
    // (cholInv ch3).to_string(),
    ch3.to_string(),
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


