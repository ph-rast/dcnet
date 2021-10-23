/*
  Perform inverse vectorization and return a correlation matrix 
  @param V array of vector to be turned int lower cholesky of a corrleation matrix
  @return matrix[Dsd] icr
 */
matrix invvec_to_corr(vector V, int nt) {
  int vnt = num_elements(V);
  //int nt = (sqrt( 8*vnt + 9 )-1 ) / 2; // compute dimension of corrmat
  // Note that we use cholesky L to compute L*L; L[1,1] = 1. Hence, vector V fills lower tri,
  // minus that first L[1,1] element
  matrix[nt,nt] L = diag_matrix(rep_vector(1.0, nt));
  real one = 1.0;
  matrix[nt,nt] out;
  matrix[nt,nt] D;
  matrix[nt,nt] R;
  vector[vnt+1] rV = append_row(one, V);
  int index = 0;
  for(i in 1:nt) {
    for(j in 1:nt) { 
      if( i >= j ) {
	index = index + 1;
	L[i,j] = rV[index];
      }
      if (i <j) {
	L[i,j] = 0.0;
	  }
    } // invvec
  } // invvec
  out =  L*L' ;
  D = diag_matrix( inv_sqrt( diagonal(out) ) );
  R = D * out * D;
  return(R);
}
