summary.dcnet <- function(object   ) {
  dt_object <- as.data.table(object$model_fit$summary( ))
  
  ## Parameters for each model
  common_params <- c("lp")
  var_params <- c("phi0_fixed|phi0_L|phi0_tau|vec_phi_fixed|phi_L|phi_tau")

  ccc_params <- c("c_h", "a_h", "b_h", "R", "beta", "c_h_var")
  dcc_D_params <- c("c_h_fixed|a_h_fixed|b_h_fixed|c_h_L|c_h_tau|a_h_L|a_h_tau|b_h_L|b_h_tau")
  dcc_Q_params <- c("l_a_q|l_a_q_r|l_b_q|l_b_q_r|S")
  dcc_params <- paste0(dcc_D_params, '|', dcc_Q_params)   
}
