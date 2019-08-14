
use crate::Number;

use crate::cov::NA;

use std::io::Write;

/// CHOLESKY DECOMPOSITION

///   Simple Cholesky decomposition of a symmetric, positive definite matrix.
///   The result for a matrix  M  is a lower triangular matrix  L  such that:
///
///      M = L L^T
///
///   Example:
///
/// >            (  2 -1  0 )   (  1.41  0     0    )
/// >            ( -1  2 -1 )   ( -0.70  1.22  0    )
/// > choldx     (  0 -1  2 ) = (  0.00 -0.81  1.15 )
///
/// Given a positive-definite symmetric matrix a[1..n][1..n],
/// this routine constructs its Cholesky decomposition,
/// A = L · L^T
/// The Cholesky factor L is returned in the lower triangle of a,
/// except for its diagonal elements which are returned in p[1..n].

pub fn do_choldc(a: &mut NA, n: usize) {

    let ll = n*n;
    let mut arr: [Number; 31] = [0.0; 31];

// -- access to arrays of symmetrical matrices
// indV :: Int -> Int -> Int -> Int
// indV w i0 j0 = i0*w+j0 -- w=nj width of niXnj matrix, i0=0..ni-1, j0=0..nj-1
// indVs :: Int -> Int -> Int -> Int
// indVs w i0 j0 | i0 <= j0   = i0*w - (i0*(i0-1)) `div` 2 + j0-i0
//               | otherwise = j0*w - (j0*(j0-1)) `div` 2 + i0-j0
// let ixa = indVs n
//     ixarr = indV n
    let w = &n;
    let ixa = |i0: usize, j0: usize| {
        if i0 <= j0 { j0 + i0*w - (i0*(i0+1))/2 }  else { i0 + j0*w - (j0*(j0+1))/2 }
    };
    let ixarr = |i0: usize, j0: usize| i0*w+j0;

    for i0 in 0..n {
        for j0 in i0..n {
            let aij = a[ixa(i0, j0)];
            if i0 == j0 { arr[ll + i0] = aij; } else { arr[ixarr(j0, i0)] = aij; };
            for k0 in 0..=i0 {
                let aik = arr[ixarr(i0, k0)];
                let ajk = arr[ixarr(j0, k0)];
                let maij = if i0 == j0 { arr[ll+i0] } else { arr[ixarr(j0, i0)] };
                let s = maij - aik*ajk;
                if i0 == j0 { arr[ll+i0]=s; } else { arr[ixarr(j0, i0)] = s; };
            };
            let msum = if i0 == j0 { arr[ll+i0] } else { arr[ixarr(j0, i0)] };
            let s = if i0==j0 && msum < 0.0 {
                writeln!(std::io::stderr(), 
                         "choldc: not a positive definite matrix ")
                    .unwrap();
                std::process::exit(1);
            } else { msum };
            let p_ip = arr[ll+i0];
            let p = if i0==j0 { s.sqrt() } else { s/p_ip };
            if i0==j0 { arr[ll+i0] = p; } else { arr[ixarr(j0, i0)] = p; };
        }
    }
    // -- copy diagonal back into array
    for i0 in 0..n {
      let aii = arr[ll+i0];
      arr[ixarr(i0, i0)] = aii;
    }
    // for i0 in 0..ll { vj[i0] = arr[i0]; };
    a.clone_from_slice(&arr[..ll]);
}

///   Matrix inversion using Cholesky decomposition of a symmetric, positive definite matrix.
///   based on Numerical Recipies formula in 2.9
///   The result for a matrix  M  is a lower triangular matrix  L  such that:
///
///   Example:
///
/// >            (  2 -1  0 )   (  0.75  0.50  0.25 )
/// >            ( -1  2 -1 )   (  0.50  1.00  0.50 )
/// > cholinv    (  0 -1  2 ) = (  0.25  0.50  0.75 )
///
pub fn do_cholinv(a: &mut NA, n: usize) {

    let ll = n*n;
    let mut arr: [Number; 31] = [0.0; 31];
    let w = &n;
    let ixa = |i0: usize, j0: usize| {
        if i0 <= j0 {
            if i0 == 0 { j0 } else { i0*w - (i0*(i0-1)) / 2 + j0-i0 }
        }  else {
            if j0 == 0 { i0 } else { j0*w - (j0*(j0-1)) / 2 + i0-j0 }
        }
    };
    let ixarr = |i0: usize, j0: usize| i0*w+j0;

    for i0 in 0..n {
        for j0 in i0..n {
            let aij = a[ixa(i0, j0)];
            if i0 == j0 { arr[ll + i0] = aij; } else { arr[ixarr(j0, i0)] = aij; };
            for k0 in 0..=i0 {
                let aik = arr[ixarr(i0, k0)];
                let ajk = arr[ixarr(j0, k0)];
                let maij = if i0 == j0 { arr[ll+i0] } else { arr[ixarr(j0, i0)] };
                let s = maij - aik*ajk;
                if i0 == j0 { arr[ll+i0]=s; } else { arr[ixarr(j0, i0)] = s; };
            };
            let msum = if i0 == j0 { arr[ll+i0] } else { arr[ixarr(j0, i0)] };
            let s = if i0==j0 && msum < 0.0 {
                writeln!(std::io::stderr(),
                         "cholinv: not a positive definite matrix ")
                    .unwrap();
                std::process::exit(1);
            } else { msum };
            let p_ip = arr[ll+i0];
            let p = if i0==j0 { s.sqrt() } else { s/p_ip };
            if i0==j0 { arr[ll+i0] = p; } else { arr[ixarr(j0, i0)] = p; };
        }
    }


    // -- copy diagonal back into array
    for i0 in 0..n {
      let p_i = arr[ll+i0];
      arr[ixarr(i0, i0)] = 1.0/p_i;
      for j0 in i0+1..n {
          arr[ll+n] = 0.0;
          for k0 in i0..j0 {
              let ajk = arr[ixarr(j0, k0)];
              let aki = arr[ixarr(k0, i0)];
              let s = arr[ll+n];
              arr[ll+n] = s - ajk * aki;
          }

        let msum = arr[ll+n];
        let p_j = arr[ll+j0];
        arr[ixarr(j0, i0)] = msum/p_j
      }
    }

    let idx = |i0: usize, j0: usize| i0*n+j0;
    for i0 in 0..n {
    for j0 in i0..n {
        let mut aij = 0.0;
        for k0 in 0..n {
            aij += arr[idx(k0, i0)] * arr[idx(k0, j0)];
        }
        a[ixa(i0, j0)] = aij;
    }}
}
// C version Numerical Recipies 2.9
// for (i=1;i<=n;i++) {
//   for (j=i;j<=n;j++) {
//     for (s=a[i][j],k=i-1;k>=1;k--) s -= a[i][k]*a[j][k];
//     if (i == j) {
//       if (s <= 0.0) nrerror("choldc failed");
//       p[i]=sqrt(s);
//     } else a[j][i]=s/p[i];
//   }
// }
//  In this, and many other applications, one often needs L^(−1) . The lower
//  triangle of this matrix can be efficiently found from the output of choldc:
// for (i=1;i<=n;i++) {
//   a[i][i]=1.0/p[i];
//   for (j=i+1;j<=n;j++) {
//     s=0.0;
//     for (k=i;k<j;k++) s -= a[j][k]*a[k][i];
//     a[j][i]=s/p[j];
//   }
// }

