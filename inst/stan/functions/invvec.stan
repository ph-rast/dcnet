/*
  Perform inverse vectorization and return a correlation matrix 
  @param V array of SD vector to be turned into lower cholesky of a covariance matrix
  @return matrix[Dsd] icr
 */
matrix invvec_to_corr(vector V, int nt) {
  int vnt = num_elements(V);
  matrix[nt,nt] L = diag_matrix(rep_vector(1.0, nt));
  matrix[nt,nt] out;
  vector[nt] D;
  matrix[nt,nt] R;
  int index = 0;
  for(i in 1:nt) {
    for(j in 1:nt) { 
      if( i >= j ) {
	index = index + 1;
	L[i,j] = V[index]; //V contains SD's, but here we create chol of Covmat
      }
      if (i <j) {
	L[i,j] = 0.0;
	  }
    } // invvec
  } 
  out = tcrossprod(L); //L*L' ;
  D = inv_sqrt( diagonal(out) );
  R = quad_form_diag(out, D); //
  return(R);
}
