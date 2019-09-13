
use crate::*;
use crate::chol::*;

use std::fmt;

use std::convert::From;

pub type NA = [Number];
pub type NA3 = [Number; 3];
pub type NA4 = [Number; 4];
pub type NA5 = [Number; 5];
pub type NA6 = [Number; 6];
pub type NA9 = [Number; 9];
pub type NA10 = [Number; 10];
pub type NA12 = [Number; 12];
pub type NA15 = [Number; 15];
pub type NA16 = [Number; 16];
pub type NA25 = [Number; 25];

enum Dim { Dim3, Dim4, Dim5 }

#[derive(Default, Debug, PartialEq, Clone)]
pub struct Vecn<T> { pub v: T }

pub type Vec3 = Vecn<NA3>;
pub type Vec4 = Vecn<NA4>;
pub type Vec5 = Vecn<NA5>;

impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vec3:{}", pretty_matrix(3,1,&self.v))
    }
}
impl fmt::Display for Vec4 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vec4:{}", pretty_matrix(4,1,&self.v))
    }
}
impl fmt::Display for Vec5 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vec5:{}", pretty_matrix(5,1,&self.v))
    }
}
impl From<NA3> for Vec3 {
    fn from(v: NA3) -> Self {
        Vec3 { v }
    }
}
impl From<NA4> for Vec4 {
    fn from(v: NA4) -> Self {
        Vec4 { v }
    }
}
impl From<NA5> for Vec5 {
    fn from(v: NA5) -> Self {
        Vec5 { v }
    }
}
impl From<Vec<Number>> for Vec3 {
    fn from(v: Vec<Number>) -> Self {
        let a: NA3 =if v.len() >= 3 { [ v[0], v[1], v[2], ] } else { [0.0;3] };
        Vec3 { v: a }
    }
}
impl From<Vec<Number>> for Vec4 {
    fn from(v: Vec<Number>) -> Self {
        let a: NA4 =if v.len() >= 4 { [ v[0], v[1], v[2], v[3], ] } else { [0.0;4] };
        Vec4 { v: a }
    }
}
impl From<Vec<Number>> for Vec5 {
    fn from(v: Vec<Number>) -> Self {
        let a: NA5 =if v.len() >= 5 { [ v[0], v[1], v[2], v[3],  v[4], ] } else { [0.0;5] };
        Vec5 { v: a }
    }
}
// pub type Vecxx = struct Vecx3 { pub v: NA3 }

// pub enum Dim { Dim3, Dim4, Dim5 }
enum WebEvent {
    // An `enum` may either be `unit-like`,
    PageLoad,
    PageUnload,
    // like tuple structs,
    KeyPress(char),
    Paste(String),
    // or like structures.
    Click { x: i64, y: i64 },
}

#[derive(Default, Debug, PartialEq, Clone)]
pub struct Cov<T> { pub v: T }
pub type Cov3 = Cov<NA6>;
pub type Cov4 = Cov<NA10>;
pub type Cov5 = Cov<NA15>;

impl Cov3 {
    // SymMat
    pub fn diag(&self) -> NA3 {
        [self.v[0], self.v[3], self.v[5], ]
    }
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
        let xx: &mut NA9 = &mut [0.0; 9];
        xx[..6].copy_from_slice(&self.v);
        do_choldc(&mut xx[..], 3);
        Jac33 { v: *xx }
    }

    pub fn cholinv(&self) -> Cov3 {
        let xx: &mut NA6 = &mut [0.0; 6];
        xx.copy_from_slice(&self.v);
        do_cholinv(&mut xx[..], 3);
        Cov3 { v: *xx }
    }
    pub fn scale_diag(&self, s: f64) -> Cov3 {
        // Cov {v: [self.v[0]*s, self.v[1], self.v[2], self.v[3]*s, self.v[4], self.v[5]*s, ]}
        Cov {v: [self.v[0]*s, 0f64, 0f64, self.v[3]*s, 0f64, self.v[5]*s, ]}
    }
}

impl fmt::Display for Cov3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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
        write!(f, "Cov3:{}", pretty_matrix(3,3,&v))
    }
}

impl From<NA6> for Cov3 {
    fn from(v: NA6) -> Self {
        Cov3 { v }
    }
}
impl From<NA9> for Cov3 {
    fn from(v: NA9) -> Self {
        Cov3 { v: [ v[0], v[1], v[2], v[4], v[5], v[8], ] }
    }
}
impl From<Vec<Number>> for Cov3 {
    fn from(v: Vec<Number>) -> Self {
        let a: NA6 =
            if v.len() == 6 { [ v[0], v[1], v[2], v[3], v[4], v[5], ] }
            else if v.len() >= 9 { [ v[0], v[1], v[2], v[4], v[5], v[8], ] }
            else { [0.0;6] }; // error
        Cov3 { v: a }
    }
}

impl Cov4 {
    pub fn diag(&self) -> NA4 {
        [self.v[0], self.v[4], self.v[7], self.v[9], ]
    }
}
impl fmt::Display for Cov4 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let w = 4;
        let idx = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let mut v: NA16 = [0.0; 16];
        for i in 0usize..4 {
            for j in 0usize..4 {
                v[j*w+i] = self.v[idx(i,j)];
            }
        }
        write!(f, "Cov4:{}", pretty_matrix(4,4,&v))
    }
}

impl From<NA10> for Cov4 {
    fn from(v: NA10) -> Self {
        Cov4 { v }
    }
}
impl From<Vec<Number>> for Cov4 {
    fn from(v: Vec<Number>) -> Self {
        let a: NA10 =
            if v.len() == 10 { [ v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], ] }
            else if v.len() >= 16 { [ v[0], v[1], v[2], v[3], v[5], v[6], v[7], v[10], v[11], v[15], ] }
            else { [0.0;10] }; // error
        Cov4 { v: a }
    }
}
impl Cov5 {
    pub fn diag(&self) -> NA5 {
        [self.v[0], self.v[5], self.v[9], self.v[12], self.v[14], ]
    }
}
impl From<&Cov5> for Cov3 { // we make a Cov3 from a Cov5 by just dropping the last R indices...
    fn from(cv: &Cov5) -> Self {
        let v = &cv.v;
        let a: NA6 = { [ v[0], v[1], v[2], v[5], v[6], v[9], ] };
        Cov3 { v: a }
    }
}
#[derive(Default, Debug, PartialEq, Clone)]
pub struct Jac<T> { pub v: T }
pub type Jac33 = Jac<NA9>;
pub type Jac34 = Jac<NA12>;
pub type Jac35 = Jac<NA15>;
pub type Jac53 = Jac<NA15>;
pub type Jac55 = Jac<NA25>;

impl Jac33 {
    pub fn tr(&self) -> Jac33 {
        let w = 3;
        let ixa = |i0: usize, j0: usize| i0*w+j0;

        let xx: &mut NA9 = &mut [0.0; 9];
        for i0 in 0..w {
            for j0 in 0..w {
                xx[ixa(i0, j0)] = self.v[ixa(j0, i0)];
            }
        }
        Jac33 { v: *xx }
    }
}
impl From<NA9> for Jac33 {
    fn from(v: NA9) -> Self {
        Jac33 { v }
    }
}

impl fmt::Display for Jac33 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Jac33:{}", pretty_matrix(3,3,&self.v))
    }
}

impl Jac55 {
    pub fn tr(&self) -> Jac55 {
        let w = 5;
        let ixa = |i0: usize, j0: usize| i0*w+j0;

        let xx: &mut NA25 = &mut [0.0; 25];
        for i0 in 0..w {
            for j0 in 0..w {
                xx[ixa(i0, j0)] = self.v[ixa(j0, i0)];
            }
        }
        Jac55 { v: *xx }
    }
}
impl fmt::Display for Jac55 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Jac55:{}", pretty_matrix(5,5,&self.v))
    }
}
impl From<NA25> for Jac55 {
    fn from(v: NA25) -> Self {
        Jac55 { v }
    }
}
impl From<Vec<Number>> for Jac55 {
    fn from(v: Vec<Number>) -> Self {
        let a: NA25 =
            if v.len() >= 25 { [ v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15], v[16], v[17], v[18], v[19],v[20], v[21], v[22], v[23], v[24], ] }
            else { [0.0;25] }; // error
        Jac55 { v: a }
    }
}

impl Jac34 {

}
impl From<NA12> for Jac34 {
    fn from(v: NA12) -> Self {
        Jac34 { v }
    }
}

impl fmt::Display for Jac34 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Jac34:{}", pretty_matrix(3,4,&self.v))
    }
}

impl Jac53 {
    // pub fn tr(&self) -> Jac35 {
    //     let ixa = |i0: usize, j0: usize| i0*3+j0;
    //     let ixb = |i0: usize, j0: usize| i0*5+j0;

    //     let xx: &mut NA15 = &mut [0.0; 15];
    //     for i0 in 0..3 {
    //         for j0 in 0..5 {
    //             xx[ixb(i0, j0)] = self.v[ixa(j0, i0)];
    //         }
    //     }
    //     Jac35 { v: *xx }
    // }
}
impl From<NA15> for Jac53 {
    fn from(v: NA15) -> Self {
        Jac53 { v }
    }
}

impl fmt::Display for Jac53 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Jac53:{}", pretty_matrix(5,3,&self.v))
    }
}

impl From<&[Number]> for Cov5 {
    fn from(v: &[Number]) -> Self {
        let a: NA15 =
            if v.len() == 15 { [ v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], ] }
            else if v.len() >= 25 { [ v[0], v[1], v[2], v[3], v[4], v[6], v[7], v[8], v[9], v[12], v[13], v[14], v[18], v[19], v[24],] }
            else { [0.0;15] }; // error
        Cov5 { v: a }
    }
}
impl From<Vec<Number>> for Cov5 {
    fn from(v: Vec<Number>) -> Self {
        let a: NA15 =
            if v.len() == 15 { [ v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], ] }
            else if v.len() >= 25 { [ v[0], v[1], v[2], v[3], v[4], v[6], v[7], v[8], v[9], v[12], v[13], v[14], v[18], v[19], v[24],] }
            else { [0.0;15] }; // error
        Cov5 { v: a }
    }
}

impl Cov5 {
    // SymMat
    pub fn det(&self) -> Number {
        // [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o] = self.0;
        let (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) = (self.v[0],self.v[1],self.v[2],self.v[3],self.v[4],self.v[5],self.v[6],self.v[7],self.v[8],self.v[9],self.v[10],self.v[11],self.v[12],self.v[13],self.v[14],);
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

    pub fn choldc(&self) -> Jac55 {
        let xx: &mut NA25 = &mut [0.0; 25];
        xx[..15].copy_from_slice(&self.v);
        do_choldc(&mut xx[..], 5);
        Jac55 { v: *xx }
    }

    pub fn cholinv(&self) -> Cov5 {
        let xx: &mut NA15 = &mut [0.0; 15];
        xx.copy_from_slice(&self.v);
        do_cholinv(&mut xx[..], 5);
        Cov5 { v: *xx }
    }
}
impl fmt::Display for Cov5 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let w = 5;
        let idx = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let mut v: NA25 = [0.0; 25 ];
        for i in 0usize..5 {
            for j in 0usize..5 {
                v[j*w+i] = self.v[idx(i,j)];
            }
        }
        write!(f, "Cov5:{}", pretty_matrix(5,5,&v))
    }
}


/// pretty print of matrix
pub fn pretty_matrix(r: usize, c: usize, v: &[Number]) -> String {
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

use std::ops::Add;
impl Add<&Vec3> for &Vec3 {
    type Output = Vec3;
    fn add(self, other: &Vec3) -> Vec3 {
        let mut r: NA3 = [0f64; 3];
        for i in 0..3 { r[i] = self.v[i] + other.v[i]; }
        Vec3 { v: r }
    }
}
impl Add<&Vec4> for &Vec4 {
    type Output = Vec4;
    fn add(self, other: &Vec4) -> Vec4 {
        let mut r: NA4 = [0f64; 4];
        for i in 0..4 { r[i] = self.v[i] + other.v[i]; }
        Vec4 { v: r }
    }
}
impl Add<&Vec5> for &Vec5 {
    type Output = Vec5;
    fn add(self, other: &Vec5) -> Vec5 {
        let mut r: NA5 = [0f64; 5];
        for i in 0..5 { r[i] = self.v[i] + other.v[i]; }
        Vec5 { v: r }
    }
}
impl Add<&Cov3> for &Cov3 {
    type Output = Cov3;
    fn add(self, other: &Cov3) -> Cov3 {
        let mut r: NA6 = [0f64; 6];
        for i in 0..6 { r[i] = self.v[i] + other.v[i]; }
        Cov3 { v: r }
    }
}
impl Add<&Cov4> for &Cov4 {
    type Output = Cov4;
    fn add(self, other: &Cov4) -> Cov4 {
        let mut r: NA10 = [0f64; 10];
        for i in 0..10 { r[i] = self.v[i] + other.v[i]; }
        Cov4 { v: r }
    }
}
//-------------------------------------------------------------------------------
use std::ops::Sub;
impl Sub<&Vec3> for &Vec3 {
    type Output = Vec3;
    fn sub(self, other: &Vec3) -> Vec3 {
        let mut r: NA3 = [0f64; 3];
        for i in 0..3 { r[i] = self.v[i] - other.v[i]; }
        Vec3 { v: r }
    }
}
impl Sub<&Vec5> for &Vec5 {
    type Output = Vec5;
    fn sub(self, other: &Vec5) -> Vec5 {
        let mut r: NA5 = [0f64; 5];
        for i in 0..5 { r[i] = self.v[i] - other.v[i]; }
        Vec5 { v: r }
    }
}
impl Sub<&Cov3> for &Cov3 {
    type Output = Cov3;
    fn sub(self, other: &Cov3) -> Cov3 {
        let mut r: NA6 = [0f64; 6];
        for i in 0..6 { r[i] = self.v[i] - other.v[i]; }
        Cov3 { v: r }
    }
}
impl Sub<&Cov5> for &Cov5 {
    type Output = Cov5;
    fn sub(self, other: &Cov5) -> Cov5 {
        let mut r: NA15 = [0f64; 15];
        for i in 0..15 { r[i] = self.v[i] - other.v[i]; }
        Cov5 { v: r }
    }
}

//-------------------------------------------------------------------------------
use std::ops::Mul;

impl Mul<&Vec3> for &Cov3 {     // Cov3 * Vec3 -> Vec3
    type Output = Vec3;
    fn mul(self, other: &Vec3) -> Vec3 {
        //self.v.len() is 9;
        let n = 3;
        let w = n;
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };

        let mut r: NA3 = [0f64; 3];
        for i in 0..n {
        r[i] = 0.0;
        for k in 0..n {
            r[i] += self.v[ixa(i,k)] * other.v[k];
        }
        }
        Vec3 { v: r }
    }
}
impl Mul<&Vec3> for &Jac53 {     // Jac53 * Vec3 -> Vec5
    type Output = Vec5;
    fn mul(self, other: &Vec3) -> Vec5 {
        let m = 5; // corrected, mxn * nx1 -> mx1; column index counts fastest
        let n = 3;
        let ixa = |i0, j0| i0*n+j0; // = indV n

        let mut r: NA5 = [0f64;  5];
        for i in 0..m {
            r[i] = 0f64;
            for k in 0..n {
                r[i] += self.v[ixa(i,k)] * other.v[k];
            }
        }
        Vec5 { v: r }
    }
}

impl Mul<&Vec5> for &Jac53 {     // Jac53T * Vec5 -> Vec3
    type Output = Vec3;
    fn mul(self, other: &Vec5) -> Vec3 {
        let m = 5; // corrected
        let n = 3;
        let ixa = |i0, j0| i0*n+j0; // = indV n

        let mut r: NA3 = [0f64; 3];
        for i in 0..n {
            r[i] = 0f64;
            for k in 0..m {
                r[i] += self.v[ixa(k,i)] * other.v[k];
            }
        }
        Vec3 { v: r }
    }
}
impl Mul<&Vec5> for &Cov5 {     // Cov5 * Vec5 -> Vec5
    type Output = Vec5;
    fn mul(self, other: &Vec5) -> Vec5 {
        let n = 5;
        let w = n;
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };

        let mut r: NA5 = [0f64; 5];
        for i in 0..n {
        r[i] = 0.0;
        for k in 0..n {
            r[i] += self.v[ixa(i,k)] * other.v[k];
        }
        }
        Vec5 { v: r }
    }
}

impl Mul<&Vec3> for &Vec3 {     // Vec3 * Vec3 -> Number
    type Output = Number;
    fn mul(self, other: &Vec3) -> Number {
        let nb = 3;
        let mut s = 0.0;
        for k in 0..nb {
            s += self.v[k] * other.v[k];
        }
        s
    }
}
impl Mul<&Vec5> for &Vec5 {     // Vec5 * Vec5 -> Number
    type Output = Number;
    fn mul(self, other: &Vec5) -> Number {
        let nb = 5;
        let mut s = 0.0;
        for k in 0..nb {
            s += self.v[k] * other.v[k];
        }
        s
    }
}
// 
impl Mul<&Jac35> for &Cov3 {    // Cov3 * Jac53T -> Jac53T
    type Output = Jac35;
    fn mul(self, other: &Jac35) -> Jac35 {
        let m = 3; // C*J -> mxm * mxn -> mxn
        let n = 5;
        let mut r = [0f64; 15];
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*m - (i0*(i0+1))/2 }  else { i0 + j0*m - (j0*(j0+1))/2 }
        };
        let ixb = |j0, i0| i0*m+j0; // = indV n
        for i in 0..m {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..m {
                    s += self.v[ixa(i,k)] * other.v[ixb(k,j)];
                }
                r[ixb(i,j)] = s;
            }
        }
        Jac35{v: r}
    }
}
impl Mul<&Cov3> for &Jac53 {    // Jac53 * Cov3 -> Jac53
    type Output = Jac53;
    fn mul(self, other: &Cov3) -> Jac53 {
        let m = 5; // J*C -> mxn * nxn -> mxn
        let n = 3;
        let mut r = [0f64; 15];
        let ixa = |i0, j0| i0*n+j0; // = indV n
        // let w = n; //= indVs n
        let ixb = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*n - (i0*(i0+1))/2 }  else { i0 + j0*n - (j0*(j0+1))/2 }
        };
        for i in 0..m {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.v[ixa(i,k)] * other.v[ixb(k,j)];
                }
                r[ixa(i,j)] = s;
            }
        }
        Jac53{v: r}
    }
}
impl Mul<&Cov5> for &Jac35 {    // Jac35 * Cov5-> Jac35
    type Output = Jac35;
    fn mul(self, other: &Cov5) -> Jac35 {
        let m = 3; // J*C -> mxn * nxn -> mxn
        let n = 5;
        let mut r = [0f64; 15];
        let ixa = |j0, i0| i0*m+j0; // = indV n
        // let w = n; //= indVs n
        let ixb = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*n - (i0*(i0+1))/2 }  else { i0 + j0*n - (j0*(j0+1))/2 }
        };
        for i in 0..m {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.v[ixa(i,k)] * other.v[ixb(k,j)];
                }
                r[ixa(i,j)] = s;
            }
        }
        Jac35{v: r}
    }
}
impl Mul<&Jac53> for &Jac35 {    // Jac35 * Jac53-> Jac33
    type Output = Jac33;
    fn mul(self, other: &Jac53) -> Jac33 {
        let m = 3_usize; // J*J -> mxn * nxm -> mxm
        let n = 5;
        let mut r = [0f64; 9];
        let ixa = |j0, i0| i0*m+j0; // = indV n
        let ixb = |i0, j0| i0*m+j0; // = indV m
        let ixr = |i0, j0| i0*m+j0; // = indV m
        for i in 0..m {
            for j in 0..m {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.v[ixa(i,k)] * other.v[ixb(k,j)];
                }
                r[ixr(i,j)] = s;
            }
        }
        // assert_eq!(r[1], r[3]);
        // assert_eq!(r[2], r[6]);
        // assert_eq!(r[5], r[7]);
        // r.into()
        Jac33 { v: r }
    }
}



// sandwich operators, these are the two-sided Mul operator, J*C -> JT.C.J,  or V*C -> VT.C.V

use std::ops::Rem;
impl Rem<&Cov3> for &Jac34 {
    type Output = Cov4;
    fn rem(self, other: &Cov3) -> Cov4 {
        // let l = other.v.len();
        let n = 3; // match l 6->3, 10->4, 15->5
        let m = self.v.len() / n; // mxn * nxn * nxm -> mxm
        let mut vint = self.v.clone(); // mxn * nxn -> mxn, v has same size as self.v
        let w = n; //= indVs n
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let ixb = |i0, j0| i0*m+j0; // = indV m
        for i in 0..n {
        for j in 0..m {
        let mut s = 0.0;
        for k in 0..n {
            s += other.v[ixa(i,k)] * self.v[ixb(k,j)];
        }
        vint[ixb(i,j)] = s;
        }}
        let mut vp = [0.0_f64; 10];  // match l 6->3, 10->4, 15->5
        let ixi = |i0, j0| i0*m+j0; // = indV m
        let w = m; //= indVs m
        let ixc = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        for i in 0..m {
        for j in i..m {
        let mut s = 0.0;
        for k in 0..n {
            s += self.v[ixi(k,i)] * vint[ixi(k,j)];
        }
        vp[ixc(i,j)] = s;
        }}
        Cov4{v: vp}
    }
}

impl Rem<&Cov5> for &Jac55 {
    type Output = Cov5;
    fn rem(self, other: &Cov5) -> Cov5 {
        // let l = other.v.len();
        let n = 5; // match l  6->3, 10->4, 15->5
        let m = self.v.len() / n; // mxn * nxn * nxm -> mxm
        let mut vint = self.v.clone(); // mxn * nxn -> mxn, v has same size as self.v
        let w = n; //= indVs n
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let ixb = |i0, j0| i0*m+j0; // = indV m
        for i in 0..n {
        for j in 0..m {
        let mut s = 0.0;
        for k in 0..n {
            s += other.v[ixa(i,k)] * self.v[ixb(k,j)];
        }
        vint[ixb(i,j)] = s;
        }}
        let mut vp = [0f64; 15];  // match l 6->3, 10->4, 15->5
        let ixi = |i0, j0| i0*m+j0; // = indV m
        let w = m; //= indVs m
        let ixc = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        for i in 0..m {
        for j in i..m {
        let mut s = 0.0;
        for k in 0..n {
            s += self.v[ixi(k,i)] * vint[ixi(k,j)];
        }
        vp[ixc(i,j)] = s;
        }}
        Cov5{v: vp}
    }
}

impl Rem<&Cov5> for &Jac53 {
    type Output = Cov3;
    fn rem(self, other: &Cov5) -> Cov3 {
        // J53(d5/d3) C5 -> C3 by J53T * C5 * J53 -> C3
        let n = 3; // corrected
        let m = 5; // JT*C*J -> nxm * mxm -> nxm * mxn -> nxn
        let mut vint = self.v.clone(); // nxm * mxm -> nxm, vint has same size as self.v
        // let w = m; //= indVs m
        let ixb = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*m - (i0*(i0+1))/2 }  else { i0 + j0*m - (j0*(j0+1))/2 }
        };
        let ixa = |i0, j0| i0*n+j0; // = indV n T
        let ixi = |i0, j0| i0*m+j0; // indV m
        for i in 0..n {
            for j in 0..m {
                let mut s = 0.0;
                for k in 0..m {
                    s += self.v[ixa(k,i)] * other.v[ixb(k,j)];
                }
                vint[ixi(i,j)] = s;
            }
        }
        let ixa = |i0, j0| i0*n+j0; // = indV n
        let ixr = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*n - (i0*(i0+1))/2 }  else { i0 + j0*n - (j0*(j0+1))/2 }
        };
        let mut r = [0f64; 6];  // match l 6->3, 10->4, 15->5
        for i in 0..n {
            for j in i..n {
                let mut s = 0.0;
                for k in 0..m {
                    s += vint[ixi(i,k)] * self.v[ixa(k,j)];
                }
                r[ixr(i,j)] = s;
            }
        }
        Cov3{v: r}
    }
}
impl Rem<&Cov3> for &Jac53 {    // Jac53.Cov3.Jac53T -> Cov5
    type Output = Cov5;
    fn rem(self, other: &Cov3) -> Cov5 {
        // C3 -> C5 by J53 * C3 * J53T -> C5
        let n = 5; // corrected
        let m = 3; // J*C*JT -> nxm * mxm -> nxm * mxn -> nxn
        let mut vint = self.v.clone(); // nxm * mxm -> nxm, vint has same size as self.v
        // let w = m; //= indVs m
        let ixb = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*m - (i0*(i0+1))/2 }  else { i0 + j0*m - (j0*(j0+1))/2 }
        };
        let ixa = |i0, j0| i0*m+j0; // = indV m
        let ixi = |i0, j0| i0*m+j0; // indV m
        for i in 0..n {
            for j in 0..m {
                let mut s = 0.0;
                for k in 0..m {
                    s += self.v[ixa(i,k)] * other.v[ixb(k,j)];
                }
                vint[ixi(i,j)] = s;
            }
        }
        let ixa = |i0, j0| i0*m+j0; // = indV m T
        let ixr = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*n - (i0*(i0+1))/2 }  else { i0 + j0*n - (j0*(j0+1))/2 }
        };
        let mut r = [0f64; 15];  // match l 6->3, 10->4, 15->5
        for i in 0..n {
            for j in i..n {
                let mut s = 0.0;
                for k in 0..m {
                    s += vint[ixi(i,k)] * self.v[ixa(j,k)];
                }
                r[ixr(i,j)] = s;
            }
        }
        Cov5{v: r}
    }
}
impl Rem<&Cov3> for &Jac33 {    // Jac33T.Cov3.Jac33 -> Cov3
    type Output = Cov3;
    fn rem(self, other: &Cov3) -> Cov3 {
        // C3 -> C3 by J33T * C3 * J33 -> C3
        let n = 3; // J*C*JT -> nxn * nxn -> nxn * nxn -> nxn
        let ixa = |i0, j0| i0*n+j0; // = indV n
        let ixb = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*n - (i0*(i0+1))/2 }  else { i0 + j0*n - (j0*(j0+1))/2 }
        };
        let mut vint = [0f64; 9];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += other.v[ixb(i,k)] * self.v[ixa(k,j)];
                }
                vint[ixa(i,j)] = s;
            }
        }
        let mut r = [0f64; 6];
        for i in 0..n {
            for j in i..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.v[ixa(k,i)] * vint[ixa(k,j)];
                }
                r[ixb(i,j)] = s;
            }
        }
        Cov3{v: r}
    }
}

// this is special: CT.C.C -> C
impl Rem<&Cov3> for &Cov3 {
    type Output = Cov3;
    fn rem(self, other: &Cov3) -> Cov3 {
        let n = 3;
        let w = n;
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let ixi = |i0: usize, j0: usize| i0*w+j0;

        // Cov3 * Cov3 -> Jac33
        let mut inter = [0.0; 9];
        for i in 0..n {
        for j in 0..n {
        let mut s = 0.0;
        for k in 0..n {
            s += self.v[ixa(k,i)] * other.v[ixa(k,j)];
        }
        inter[ixi(i,j)] = s;
        }}

        let mut res = [0.0; 6];
        for i in 0..n {
        for j in 0..n {
        let mut s = 0.0;
        for k in 0..n {
            s += inter[ixi(i,k)] * self.v[ixa(k,j)];
        }
        res[ixa(i,j)] = s;
        }}
        Cov3 { v: res }
    }
}
impl Rem<&Cov5> for &Cov5 {
    type Output = Cov5;
    fn rem(self, other: &Cov5) -> Cov5 {
        let n = 5;
        let w = n;
        let ixa = |i0: usize, j0: usize| {
            if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
        };
        let ixi = |i0: usize, j0: usize| i0*w+j0;

        // Cov5 * Cov5 -> Jac55
        let mut inter = [0.0; 25];
        for i in 0..n {
        for j in 0..n {
        let mut s = 0.0;
        for k in 0..n {
            s += self.v[ixa(k,i)] * other.v[ixa(k,j)];
        }
        inter[ixi(i,j)] = s;
        }}

        let mut res = [0.0; 15];
        for i in 0..n {
        for j in 0..n {
        let mut s = 0.0;
        for k in 0..n {
            s += inter[ixi(i,k)] * self.v[ixa(k,j)];
        }
        res[ixa(i,j)] = s;
        }}
        Cov5 { v: res }
    }
}
#[test]
fn test_cov() {
    let _ch3 : Cov3 = [2.0, -1.0, 0.0, 2.0, -1.0, 2.0].into();
    // let ch3 = Cov3::from([2.0, -1.0, 0.0, 2.0, -1.0, 2.0]);
    // let ch5 = Cov5{ v: [2.0, -1.0, 0.0, 0.0, 0.0, 2.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0] };

    // let res = format!(r"
// chol: -----------------
// A = L * L^T             {}
// L                       {}
// L * L^T                 {}
// A^(-1) = L' * L'^T      {}
// A * A^(-1)              {}
// A = L * L^T             {}
// L                       {}
// L * L^T                 {}
// A^(-1) = L' * L'^T      {}
// A * A^(-1)              {}
// det this                {}
// ",
    // ch3,
    // ch3.choldc(),
    // ch3.choldc() * ch3.choldc().tr(),
    // ch3.cholinv(),
    // ch3.clone() * ch3.cholinv(),
    // ch5,
    // ch5.choldc(),
    // ch5.choldc() * ch5.choldc().tr(),
    // ch5.cholinv(),
    // &ch5 * &ch5.cholinv(),
    // ch5.det(),
    // );
    // print!("{}", res);

    // let v3 = Vec3 { v: [10.0,11.0,12.0] };
    // println!("Vec * Vec = {}\n", &v3 * &v3);
    // println!("Cov *. Cov = {}\n", (Cov3 {v: [1.0,2.0,3.0,4.0,5.0,6.0]}) * (Cov3 {v: [0.0,0.0,1.0,1.0,0.0,0.0]}));

    assert!(true, "test failed with");
}


