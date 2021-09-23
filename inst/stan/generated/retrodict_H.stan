// retrodict given distribution type
if ( distribution == 0 ){
  for (j in 1:J){
    for (t in 1:T) {
      //    rts_out[,t] = multi_normal_rng(mu[t,j,], H[t,]);
      corH[t,] = cov2cor(H[t,]);
      log_lik[t] = multi_normal_lpdf(rts[t,j,] | mu[t,], H[t,]);
    }
  }
 } else if ( distribution == 1 ) {
  for (j in 1:J) {
    for (t in 1:T) {
      //   rts_out[,t] = multi_student_t_rng(nu, mu[t,j,], H[t,]);
      corH[t,] = cov2cor(H[t,]);
      log_lik[t] = multi_student_t_lpdf(rts[t,j,] | nu, mu[t,j,], H[t,]);
    }
  }
 }
