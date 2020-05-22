// Space time 
#include <TMB.hpp>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_SCALAR(scaler);
  // Indices
  DATA_INTEGER( n_i );         // Total number of observations
  
  // Data
  DATA_VECTOR( c_i );       	// Count data
  DATA_IVECTOR( s_i );        // Station for each sample
  DATA_IVECTOR( t_i );        // Time for each sample
  DATA_IVECTOR( lat_i );
  DATA_VECTOR( lat_dist );
  DATA_MATRIX( X_xp );		            // Covariate design matrix
  DATA_INTEGER(use_p_y);
  DATA_INTEGER(use_c_y);
  DATA_INTEGER(use_p_s);
  DATA_INTEGER(use_c_s);
  DATA_INTEGER(use_p_lat);
  DATA_INTEGER(use_c_lat);
  DATA_INTEGER(use_p_yl);
  DATA_INTEGER(use_c_yl);
  
  // Fixed effects
  PARAMETER_VECTOR(beta_c);   //Station effect for counts
  PARAMETER_VECTOR(beta_p);   //Station effect for counts
  PARAMETER_VECTOR(eps_c_s);   //Station effect for counts
  PARAMETER_VECTOR(eps_p_s);   //Station effect for p
  PARAMETER_VECTOR(eps_c_y);   //Station effect for p
  PARAMETER_VECTOR(eps_p_y);   //Station effect for p
  PARAMETER_VECTOR(eps_p_lat);   //Station effect for p
  PARAMETER_VECTOR(eps_c_lat);   //Station effect for p
  PARAMETER_ARRAY(eps_p_yl);   //Station effect for p 
  PARAMETER_ARRAY(eps_c_yl);   //Station effect for p
  PARAMETER(fc_y_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_y_sd);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_y_rho);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_y_rho);             //Station s.d. for station effect for zero catches
  
  PARAMETER(fc_yl_rho_y);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_yl_rho_y);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_yl_rho_lat);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_yl_rho_lat);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_yl_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_yl_sd);             //Station s.d. for station effect for zero catches
  
  PARAMETER(fc_lat_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_lat_sd);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_lat_rho);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_lat_rho);             //Station s.d. for station effect for zero catches
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

  Type c_y_rho = f(fc_y_rho);//1/(1+exp(-1.0*fc_y_rho));
  Type p_y_rho = f(fp_y_rho);//1/(1+exp(-1.0*fp_y_rho));
  Type c_lat_rho = 1/(1+exp(-1.0*fc_lat_rho));//atan(fc_lat_rho)*2./3.154;//f(fc_lat_rho);//
  Type p_lat_rho = 1/(1+exp(-1.0*fp_lat_rho));//atan(fp_lat_rho)*2./3.154;//f(fp_lat_rho);//
  Type c_yl_rho_lat = 1/(1+exp(-1.0*fc_yl_rho_lat));//atan(fc_yl_rho_lat)*2./3.154;//f(fc_yl_rho_lat);//
  Type p_yl_rho_lat = 1/(1+exp(-1.0*fp_lat_rho));//f(fp_yl_rho_lat);//atan(fp_yl_rho_lat)*2./3.154;//
  Type c_yl_rho_y = f(fc_yl_rho_y);//1/(1+exp(-1.0*fc_yl_rho_y));//atan(fc_yl_rho_y)*2./3.154;//
  Type p_yl_rho_y = f(fp_yl_rho_y);//atan(fp_yl_rho_y)*2./3.154;//1/(1+exp(-1.0*fp_lat_rho));
  
  if(use_c_y){
    jnll_comp(0) += AR1(c_y_rho)(eps_c_y);
    jnll_comp(0) -= dnorm(eps_c_y(0),Type(0.),Type(1.),true);
  }
  if(use_p_y){
    jnll_comp(0) += AR1(p_y_rho)(eps_p_y);
    jnll_comp(0) -= dnorm(eps_p_y(0),Type(0.),Type(1.),true);
  }
  
  int n = lat_dist.size();
  matrix<Type> Sigma_c(n,n); 
  matrix<Type> Sigma_p(n,n); 
  matrix<Type> Sigma_c_yl(n,n); 
  matrix<Type> Sigma_p_yl(n,n); 
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(use_c_lat){
        Sigma_c(i,j) = pow(c_lat_rho,abs(lat_dist(i)-lat_dist(j))/scaler);
      }
      if(use_p_lat){
        Sigma_p(i,j) = pow(p_lat_rho,abs(lat_dist(i)-lat_dist(j))/scaler);
      }
      if(use_c_yl){
        Sigma_c_yl(i,j) = pow(c_yl_rho_lat,abs(lat_dist(i)-lat_dist(j))/scaler);
      }
      if(use_p_yl){
        Sigma_p_yl(i,j) = pow(p_yl_rho_lat,abs(lat_dist(i)-lat_dist(j))/scaler);
      }
    }
  }

  MVNORM_t<Type> c_dmnorm2(Sigma_c);
  MVNORM_t<Type> p_dmnorm2(Sigma_p);
  if(use_c_lat){
    jnll_comp(0) += c_dmnorm2(eps_c_lat);
  }
  if(use_p_lat){
    jnll_comp(0) += p_dmnorm2(eps_p_lat);
  }

  MVNORM_t<Type> c_dmnorm(Sigma_c_yl);
  MVNORM_t<Type> p_dmnorm(Sigma_p_yl);
  if(use_c_yl){
    jnll_comp(0) += SEPARABLE(AR1(c_yl_rho_y),c_dmnorm)(eps_c_yl);
  }
  
  if(use_p_yl){
    jnll_comp(0) += SEPARABLE(AR1(p_yl_rho_y),p_dmnorm)(eps_p_yl);
  }
  
  for (int i=0; i<n_i; i++){
    p_i(i) = 1.0/(1.0+exp(-1.0*(eta_p(i)+
      eps_p_s(s_i(i))+
      eps_p_y(t_i(i)) * exp(fp_y_sd) +
      eps_p_lat(lat_i(i)) * exp(fp_lat_sd))) +
      eps_p_yl(lat_i(i),t_i(i)) * exp(fp_yl_sd));
    
    if(c_i(i)==0.){
      jnll_comp(2) -= log(1.-p_i(i));
    } 
    
    if(c_i(i)!=0.){
      log_chat_i(i) = eta_c(i) +
        eps_c_s(s_i(i)) +
        eps_c_y(t_i(i)) * exp(fc_y_sd) +
        eps_c_lat(lat_i(i)) * exp(fc_lat_sd) +
        eps_c_yl(lat_i(i),t_i(i)) * exp(fc_yl_sd);
      
      jnll_comp(2) -= dnorm(log(c_i(i)), log_chat_i(i), exp(fsd), true );
      jnll_comp(2) -= log(p_i(i));
    }
  }
  
  jnll_comp(3) = 0.;
  for(int i=0;i<eps_c_s.size();i++){
    if(use_c_s){
      jnll_comp(3) -= dnorm(eps_c_s(i),Type(0.),exp(fs_c_sd),true);
    }
    if(use_p_s){
      jnll_comp(3) -= dnorm(eps_p_s(i),Type(0.),exp(fs_p_sd),true);
    }
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
  REPORT( Sigma_c );
  REPORT( Sigma_c );
  REPORT( Sigma_c_yl );
  REPORT( Sigma_c_yl );
  REPORT( eps_c_lat );
  REPORT( eps_p_lat );
  REPORT( eps_p_yl );
  REPORT( eps_c_yl );
  
  
  REPORT(p_i);
  
  REPORT( jnll );
  return jnll;
}
