array[J] matrix[T, nt] rts_out;
array[J] vector[T] log_lik;

if ( distribution == 0 ){
  for( j in 1:J ) {
    for (t in 1:T) {
      rts_out[j,t] = multi_normal_rng(mu[j,t], H[j,t])';
      log_lik[j,t] = multi_normal_lpdf(rts[j,t] | mu[j,t], H[j,t]);
    }
  }
 } else if ( distribution == 1 ) {
  for( j in 1:J ) {
    for (t in 1:T) {
      rts_out[j,t] = multi_student_t_rng(nu, mu[j,t], H[j,t])';
      log_lik[j,t] = multi_student_t_lpdf(rts[j,t] | nu, mu[j,t], H[j,t]);
    }
    }
 }
