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
  DATA_IVECTOR( lo_i );
  DATA_VECTOR( lat_dist );
  DATA_MATRIX( X_xp );		            // Covariate design matrix
  DATA_INTEGER(use_p_y);
  DATA_INTEGER(use_c_y);
  DATA_INTEGER(use_p_s);
  DATA_INTEGER(use_c_s);
  DATA_INTEGER(use_p_lat);
  DATA_INTEGER(use_c_lat);
  DATA_INTEGER(use_p_lo);
  DATA_INTEGER(use_c_lo);
  DATA_INTEGER(use_p_yll);
  DATA_INTEGER(use_c_yll);
  
  // Fixed effects
  PARAMETER_VECTOR(beta_c);   //Station effect for counts
  PARAMETER_VECTOR(beta_p);   //Station effect for counts
  PARAMETER_VECTOR(eps_c_s);   //Station effect for counts
  PARAMETER_VECTOR(eps_p_s);   //Station effect for p
  PARAMETER_VECTOR(eps_c_y);   //Station effect for p
  PARAMETER_VECTOR(eps_p_y);   //Station effect for p
  PARAMETER_VECTOR(eps_p_lat);   //Station effect for p
  PARAMETER_VECTOR(eps_c_lat);   //Station effect for p
  PARAMETER_VECTOR(eps_p_lo);   //Station effect for p
  PARAMETER_VECTOR(eps_c_lo);   //Station effect for p
  PARAMETER_ARRAY(eps_p_yll);   //Station effect for p 
  PARAMETER_ARRAY(eps_c_yll);   //Station effect for p
  PARAMETER(fc_y_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_y_sd);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_y_rho);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_y_rho);             //Station s.d. for station effect for zero catches
  
  PARAMETER(fc_yll_rho_y);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_yll_rho_y);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_yll_rho_lat);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_yll_rho_lat);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_yll_rho_lo);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_yll_rho_lo);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_yll_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_yll_sd);             //Station s.d. for station effect for zero catches
  
  PARAMETER(fc_lat_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_lat_sd);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_lo_sd);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_lo_sd);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_lat_rho);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_lat_rho);             //Station s.d. for station effect for zero catches
  PARAMETER(fc_lo_rho);             //Station s.d. for station effect for positive catches
  PARAMETER(fp_lo_rho);             //Station s.d. for station effect for zero catches
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
  Type c_lo_rho = 1/(1+exp(-1.0*fc_lo_rho));//atan(fc_lat_rho)*2./3.154;//f(fc_lat_rho);//
  Type p_lo_rho = 1/(1+exp(-1.0*fp_lo_rho));//atan(fp_lat_rho)*2./3.154;//f(fp_lat_rho);//
  Type c_yll_rho_lat = 1/(1+exp(-1.0*fc_yll_rho_lat));//atan(fc_yl_rho_lat)*2./3.154;//f(fc_yl_rho_lat);//
  Type p_yll_rho_lat = 1/(1+exp(-1.0*fp_yll_rho_lat));//f(fp_yl_rho_lat);//atan(fp_yl_rho_lat)*2./3.154;//
  Type c_yll_rho_lo = 1/(1+exp(-1.0*fc_yll_rho_lo));//atan(fc_yl_rho_lat)*2./3.154;//f(fc_yl_rho_lat);//
  Type p_yll_rho_lo = 1/(1+exp(-1.0*fp_yll_rho_lo));//f(fp_yl_rho_lat);//atan(fp_yl_rho_lat)*2./3.154;//
  Type c_yll_rho_y = f(fc_yll_rho_y);//1/(1+exp(-1.0*fc_yl_rho_y));//atan(fc_yl_rho_y)*2./3.154;//
  Type p_yll_rho_y = f(fp_yll_rho_y);//atan(fp_yl_rho_y)*2./3.154;//1/(1+exp(-1.0*fp_lat_rho));
  
  array<Type> pred_c_ef(eps_c_lo.size(),eps_c_lat.size(),eps_c_s.size());
  array<Type> pred_p_ef(eps_p_lo.size(),eps_p_lat.size(),eps_p_s.size());
  
  if(use_c_y){
    jnll_comp(0) += AR1(c_y_rho)(eps_c_y);
    jnll_comp(0) -= dnorm(eps_c_y(0),Type(0.),Type(1.),true);
  }
  if(use_p_y){
    jnll_comp(0) += AR1(p_y_rho)(eps_p_y);
    jnll_comp(0) -= dnorm(eps_p_y(0),Type(0.),Type(1.),true);
  }
  
  if(use_c_lat){
    jnll_comp(0) += AR1(c_lat_rho)(eps_c_lat);
  }
  if(use_p_lat){
    jnll_comp(0) += AR1(p_lat_rho)(eps_p_lat);
  }

  if(use_c_lo){
    jnll_comp(0) += AR1(c_lo_rho)(eps_c_lo);
  }
  if(use_p_lo){
    jnll_comp(0) += AR1(p_lo_rho)(eps_p_lo);
  }
  
  if(use_c_yll){
    jnll_comp(0) += AR1(c_yll_rho_y,AR1(c_yll_rho_lat,AR1(c_yll_rho_lo)))(eps_c_yll);
  }
  
  if(use_p_yll){
    jnll_comp(0) += AR1(p_yll_rho_y,AR1(p_yll_rho_lat,AR1(p_yll_rho_lo)))(eps_p_yll);
  }
  
  for (int i=0; i<n_i; i++){
    p_i(i) = 1.0/(1.0+exp(-1.0*(eta_p(i)+
      eps_p_s(s_i(i))+
      eps_p_y(t_i(i)) * exp(fp_y_sd) +
      eps_p_lat(lat_i(i)) * exp(fp_lat_sd) +
      eps_p_lo(lo_i(i)) * exp(fp_lo_sd) +
      eps_p_yll(lo_i(i),lat_i(i),t_i(i)) * exp(fp_yll_sd))));
    
    if(c_i(i)==0.){
      jnll_comp(2) -= log(1.-p_i(i));
    } 
    
    if(c_i(i)!=0.){
      log_chat_i(i) = eta_c(i) +
        eps_c_s(s_i(i)) +
        eps_c_y(t_i(i)) * exp(fc_y_sd)+
        eps_c_lat(lat_i(i)) * exp(fc_lat_sd) +
        eps_c_lo(lo_i(i)) * exp(fc_lo_sd) +
        eps_c_yll(lo_i(i),lat_i(i),t_i(i)) * exp(fc_yll_sd);
        
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

  for(int lo = 0; lo<eps_p_lo.size(); lo++){
    for(int lat = 0; lat<eps_p_lat.size(); lat++){
      for(int yy = 0; yy<eps_p_y.size(); yy++){
        pred_c_ef(lo,lat,yy) =  beta_c(0) +
          eps_c_y(yy) * exp(fc_y_sd)+
          eps_c_lat(lat) * exp(fc_lat_sd) +
          eps_c_lo(lo) * exp(fc_lo_sd) +
          eps_c_yll(lo,lat,yy) * exp(fc_yll_sd);

        pred_p_ef(lo,lat,yy) =  beta_p(0) +
          eps_p_y(yy) * exp(fp_y_sd)+
          eps_p_lat(lat) * exp(fp_lat_sd) +
          eps_p_lo(lo) * exp(fp_lo_sd) +
          eps_p_yll(lo,lat,yy) * exp(fp_yll_sd);
        
      }
    }
  }
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
  // REPORT( Sigma_c );
  // REPORT( Sigma_p );
  // REPORT( Sigma_c_yl );
  // REPORT( Sigma_p_yl );
  REPORT( eps_c_lat );
  REPORT( eps_p_lat );
  REPORT( eps_c_lo );
  REPORT( eps_p_lo );
  REPORT( eps_p_yll );
  REPORT( eps_c_yll );
  
  
  REPORT(p_i);
  REPORT(pred_c_ef);
  REPORT(pred_p_ef);
  
  REPORT( jnll );
  
  ADREPORT(eps_c_y);
  ADREPORT(eps_p_y);
  ADREPORT(eps_c_lat);
  ADREPORT(eps_p_lat);
  
  return jnll;
}
