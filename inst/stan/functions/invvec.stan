/*
  Perform inverse vectorization and return a correlation matrix via cholesky
  @param V array of SD vector to be turned into lower cholesky of a covariance matrix
  @return matrix[Dsd] icr
 */
matrix invvec_chol_to_corr(vector V, int nt) {
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
  D = inv_sqrt( diagonal(out) ); // This makes S_Lv basically a chol of cov
  R = quad_form_diag(out, D); //
  return(R);
}

/* from mc-stan discussion w. Stephen
 */
matrix convex_combine_cholesky(matrix global_chol_cor, matrix local_chol_cor, real alpha){
  int dim = rows(local_chol_cor);
  int global_dim = rows(global_chol_cor);
  matrix[global_dim,global_dim] global_cor = multiply_lower_tri_self_transpose(global_chol_cor);
  matrix[dim,dim] local_cor = multiply_lower_tri_self_transpose(local_chol_cor);
  matrix[dim,dim] L;
  matrix[dim,dim] combined_chol_cor;
  
  L = cholesky_decompose((1 - alpha)*global_cor + (alpha)*local_cor);
  combined_chol_cor=tcrossprod(L); // LL'
  return(combined_chol_cor);
}
