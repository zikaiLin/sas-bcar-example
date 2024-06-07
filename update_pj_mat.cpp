#include <RcppArmadillo.h>
#include <RcppDist.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix update_pj_mat(IntegerMatrix delta_mat, IntegerMatrix partition, IntegerVector radius_partition){
  int nr = delta_mat.nrow();
  int nc = delta_mat.ncol();
  
  NumericMatrix pj_mat(nr, nc);
  
  Rcpp::Function expGrid("expand.grid");
  
  for(int i=0; i<nr; i++){
    for(int j=0;j<nc; j++){
      int radius = radius_partition(partition(i,j)-1);
      List ij_circle_list = expGrid(seq(i-radius, i+radius), seq(j-radius, j+radius));
      IntegerVector col1 = ij_circle_list[0];
      IntegerVector col2 = ij_circle_list[1];
      IntegerMatrix ij_circle= cbind(col1, col2);
      // Rcout << ij_circle << "\n";
      int count = 0;
      double sum = 0;
      for(int s=0; s<ij_circle.nrow();s++){
        IntegerVector s_idx = ij_circle.row(s);
        if(s_idx(0) < 0 | s_idx(1) < 0 | s_idx(0) >= nr |  s_idx(1) >= nc){
          continue;
        }else{
          sum += delta_mat(s_idx(0), s_idx(1));
          count += 1;
        }
        
        pj_mat(i,j) = sum/count;
        
      }
      
    }
    
  }
  
  return pj_mat;
}


