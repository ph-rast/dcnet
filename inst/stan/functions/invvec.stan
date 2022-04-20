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
      if (i < j) {
	L[i,j] = 0.0;
	  }
    } // invvec
  } 
  out = multiply_lower_tri_self_transpose(L); //L*L' ;
  D = inv_sqrt( diagonal(out) ); // This makes S_Lv basically a chol of cov
  R = quad_form_diag(out, D); //
  return(R);
}

/* Convex combination of corr mats
   but see: https://gmarti.gitlab.io/ml/2021/02/13/swelling-effect-spd-covariance.html
   This solution leads to swelling, but should be not a problem in the current application.
 */
matrix convex_combine_cholesky(matrix global_chol_cor, matrix local_chol_cor, real alpha){
  int dim = rows(local_chol_cor);
  matrix[dim,dim] global_cor = multiply_lower_tri_self_transpose(global_chol_cor);
  matrix[dim,dim] local_cor = multiply_lower_tri_self_transpose(local_chol_cor);  
  matrix[dim,dim] combined_cor;
  
  combined_cor = (1 - alpha)*global_cor + (alpha)*local_cor;
  return(combined_cor);
}

// Swelling should not be a problem (in terms of interpretatino), as I'm assuming that
// there's a sampling distribution of correlation matrices. The average correlation
// would then reflect the fixed effet and individual deviations are random effects.
//
// Riemannian Log-Cholesky Metric: Average over log_chol, should not lead to swelling
// cf. Lin (2019) Riemannian geometry of Symmetric PD Matrices via Chol decomp.
matrix geodesic_log_chol(matrix global_chol_cor, matrix local_chol_cor, real alpha){
  int dim = rows(local_chol_cor);
  matrix[dim,dim] out;
  matrix[dim,dim] tri_global;
  matrix[dim,dim] tri_local;
  matrix[dim,dim] L;
  vector[dim] D;
  matrix[dim,dim] R;

  // Create lower tri matrix with 0 on diagonal:
  // cf. geodesic.chol function from
  // 
  tri_global=global_chol_cor - diag_matrix(diagonal(global_chol_cor));
  tri_local=local_chol_cor - diag_matrix(diagonal(local_chol_cor));;
  
  L=tri_global + alpha*(tri_local-tri_global) +
    diag_matrix(
		diagonal(global_chol_cor) .* exp(alpha * ( log(diagonal(local_chol_cor)) -
							   log(diagonal(global_chol_cor))))
		);
  
  out=multiply_lower_tri_self_transpose(L);
  D = inv_sqrt( diagonal(out) ); // Ensure matrix is a corrmat
  R = quad_form_diag(out, D); // 
  return(R);
}

matrix convex_combine(matrix global_cor, matrix local_cor, real alpha){
  int dim = rows(local_cor);
  matrix[dim,dim] R;
  matrix[dim,dim] combined_cor;
  vector[dim] D;

  combined_cor = exp( (1 - alpha)*log(global_cor) + (alpha)*log(local_cor));
  D = inv_sqrt( diagonal(combined_cor) ); // Ensure matrix is a corrmat
  R = quad_form_diag(combined_cor, D); // 
  return(R);
}
