// Space time 
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Indices
  DATA_INTEGER( n_i );         // Total number of observations
  
  // Data
  DATA_VECTOR( c_i );       	// Count data
  DATA_IVECTOR( s_i );        // Station for each sample
  DATA_IVECTOR( t_i );        // Time for each sample
  DATA_MATRIX( X_xp );		            // Covariate design matrix

  // Fixed effects
  PARAMETER_VECTOR(beta_c);   //Station effect for counts
  PARAMETER_VECTOR(beta_p);   //Station effect for counts
  PARAMETER_VECTOR(eps_c_s);   //Station effect for counts
  PARAMETER_VECTOR(eps_p_s);   //Station effect for p
  PARAMETER_VECTOR(eps_c_y);   //Station effect for p
  PARAMETER_VECTOR(eps_p_y);   //Station effect for p
  PARAMETER(fc_y_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_y_sd);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_y_rho);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_y_rho);             //Station s.d. for station effect for zero catches
  PARAMETER(fs_c_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fs_p_sd);             //Station s.d. for station effect for zero catches
  PARAMETER(fsd);
  
  // objective function -- joint negative log-likelihood
  using namespace density;
  Type jnll = 0;
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();


  //Covariate effects
  vector<Type> eta_c = X_xp * beta_c;
  vector<Type> eta_p = X_xp * beta_p;
  vector<Type> p_i(n_i);
  // Likelihood contribution from observations
  vector<Type> log_chat_i(n_i); //predicted catch
  Type p; //predicted probability of zero
  
  Type c_y_rho = 1/(1+exp(-1.0*fc_y_rho));
  Type p_y_rho = 1/(1+exp(-1.0*fp_y_rho));
  
  jnll_comp(0) += AR1(c_y_rho)(eps_c_y);
  jnll_comp(0) += AR1(p_y_rho)(eps_p_y);
  jnll_comp(0) -= dnorm(eps_p_y(0),Type(0.),Type(1.),true);
  jnll_comp(0) -= dnorm(eps_c_y(0),Type(0.),Type(1.),true);
  
  for (int i=0; i<n_i; i++){
    p_i(i) = 1.0/(1.0+exp(-1.0*(eta_p(i)+
      eps_p_s(s_i(i))+
      eps_p_y(t_i(i)) * exp(fp_y_sd))));
    
    if(c_i(i)==0.){
      jnll_comp(2) -= log(1.-p_i(i));
    } 
    
    if(c_i(i)!=0.){
      log_chat_i(i) = eta_c(i) +
        eps_c_s(s_i(i)) +
        eps_c_y(t_i(i)) * exp(fc_y_sd);
      
      jnll_comp(2) -= dnorm(log(c_i(i)), log_chat_i(i), exp(fsd), true );
      jnll_comp(2) -= log(p_i(i));
    }
  }
  
  jnll_comp(3) = 0.;
  for(int i=0;i<eps_c_s.size();i++){
    jnll_comp(3) -= dnorm(eps_c_s(i),Type(0.),exp(fs_c_sd),true);
    jnll_comp(3) -= dnorm(eps_p_s(i),Type(0.),exp(fs_p_sd),true);
  }
  
  jnll = jnll_comp.sum();

  // // Diagnostics
  REPORT( jnll_comp );
  REPORT( log_chat_i );
  REPORT( eta_c );
  REPORT( eta_p );
  REPORT( beta_c );
  REPORT( beta_p );
  REPORT( eps_p_s );
  REPORT( eps_c_s );
  REPORT( eps_c_y );
  REPORT( eps_p_y );
  
  REPORT(p_i);
  ADREPORT(eps_c_y);
  ADREPORT(eps_p_y);
  
  REPORT( jnll );
  return jnll;
}
