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
  DATA_VECTOR( c_i );       	// Count data
  DATA_IVECTOR( s_i );       	// Station
  DATA_IVECTOR( s_i_re );       	// Stations with positive counts
  DATA_IVECTOR( y_i );       	// Stations with positive counts
  DATA_VECTOR(transDist);
  DATA_IVECTOR(trans_i);
  DATA_INTEGER(est_tr);
  DATA_INTEGER(est_y);
  DATA_INTEGER(est_s);
  DATA_INTEGER(est_yt);
  
  // Fixed effects
  PARAMETER(fp);   //proportin with catches
  PARAMETER(ln_lam);
  PARAMETER(ln_theta);
  PARAMETER_VECTOR(eps_s);
  PARAMETER(lnsig_eps_s);
  PARAMETER_VECTOR(eps_y);
  PARAMETER(frho_y);
  PARAMETER(lnsig_eps_y);
  PARAMETER_VECTOR(eps_tr);
  PARAMETER(frho_tr);
  PARAMETER(lnsig_eps_tr);

  PARAMETER_ARRAY(eps_yt);
  PARAMETER(frho_yt_1);
  PARAMETER(frho_yt_2);
  PARAMETER(lnsig_eps_yt);
  
  // objective function -- joint negative log-likelihood
  using namespace density;
  Type jnll = 0;

  int n_i = c_i.size();
  
  // Likelihood contribution from observations
  Type p = exp(fp)/(1+exp(fp));
  Type rho_y = 1/(1+exp(-frho_y)); //atan(frho_y)*3.154/2.; //
  Type rho_tr = exp(frho_tr);//1/(1+exp(-frho_tr)); //atan(frho_tr)*3.154/2.;//Type(2)/(Type(1) + exp(-Type(2) * frho_y)) - Type(1);//
  Type rho_yt_1 = atan(frho_yt_1)*3.154/2.;//1/(1+exp(-frho_yt_1)); //Type(2)/(Type(1) + exp(-Type(2) * frho_yt_1)) - Type(1);//
  Type rho_yt_2 = exp(frho_yt_2);//1/(1+exp(-frho_yt_2));//atan(frho_yt_2)*3.154/2.;//
  Type theta = exp(ln_theta);
  Type sig_eps_s = exp(lnsig_eps_s);
  Type sig_eps_y = exp(lnsig_eps_y);
  Type sig_eps_tr = exp(lnsig_eps_tr);
  Type sig_eps_yt = exp(lnsig_eps_yt);
  Type a;
  
  vector<Type> log_chat_i(n_i);
  int nt = transDist.size();
  
  //Temporal 1-D
  if(est_y){
    jnll += AR1(rho_y)(eps_y);
    jnll -= dnorm(eps_y(0),Type(0.),Type(1.0),true);
  }
  
  matrix<Type> Sigma_yt(nt,nt);
  if(est_yt){
    //Spatial 1-D
    for(int i=0;i<nt;i++){
      for(int j=0;j<nt;j++){
        if(j!=i){
          Sigma_yt(i,j) = exp(-1.*rho_yt_2*abs(transDist(i)-transDist(j)));
        }
        if(j==i){
          Sigma_yt(i,j) = 1.;
        }
      }
    }
    // jnll += AR1(rho_yt_2,AR1(rho_yt_1))(eps_yt);
    jnll += SEPARABLE(AR1(rho_yt_1),MVNORM(Sigma_yt))(eps_yt);
    
    // for(int i=0;i<eps_tr.size();i++){
    //   jnll -= dnorm(eps_yt(i,0),Type(0.), Type(1.0),true);
    // }
  }
  if(est_s){
    jnll -= dnorm(eps_s,Type(0.),sig_eps_s,true).sum();
  }
  
  //Spatial 1-D
  matrix<Type> Sigma(nt,nt);
  for(int i=0;i<nt;i++){
    for(int j=0;j<nt;j++){
      if(j!=i){
        Sigma(i,j) = exp(-1.*rho_tr*abs(transDist(i)-transDist(j)));
      }
      if(j==i){
        Sigma(i,j) = 1.;
      }
    }
  }
      
  if(est_tr){
    // jnll += AR1(rho_tr)(eps_tr);
    jnll += MVNORM(Sigma)(eps_tr);
    // for(int i=1;i<nt;i++){
    //   jnll -= dnorm(eps_tr(i),pow(rho_tr,abs(transDist(i)-transDist(i-1))/0.1)*eps_tr(i-1),1.,true);
    // }
  }

    
  
  //Observation model
  for (int i=0; i<n_i; i++){
    if(c_i(i)==0.){
      jnll -= log(1.-p);
    } 
    if(c_i(i)!=0.){
      log_chat_i(i) = ln_lam + 
        eps_s(s_i_re(s_i(i))) + //Station effect
        eps_y(y_i(i)) * sig_eps_y +//Year effect
        eps_tr(trans_i(i)) * sig_eps_tr + //Transect/latitude effect
        eps_yt(trans_i(i),y_i(i))*sig_eps_yt; //Interaction between transect and latitude
      
      jnll -= log(p);
      a = pow(theta,-2);
      jnll -= dgamma(c_i(i), a, exp(log_chat_i(i))/a, true );
      
    }
  }
  
  SIMULATE {
    if(est_y){
      AR1(rho_y).simulate(eps_y);
    }else{
      eps_y.fill(0.);
    }
    if(est_tr){
      MVNORM(Sigma).simulate(eps_tr);
    }else{
      eps_tr.fill(0.);
    }
    if(est_yt){
      AR1(rho_yt_1,MVNORM(Sigma_yt)).simulate(eps_yt);
      // AR1(rho_yt_2,AR1(rho_yt_1)).simulate(eps_yt);
    }else{
      eps_yt.fill(0.);
    }
    if(est_s){
      vector<Type> eps_s_0(eps_s.size());
      eps_s_0.fill(0.);
      eps_s = rnorm(eps_s_0,sig_eps_s);
    }else{
      eps_s.fill(0.);
    }
    
    for (int i=0; i<n_i; i++){
      if(c_i(i)>0.){
        log_chat_i(i) = ln_lam + 
        eps_s(s_i_re(s_i(i))) + 
        eps_y(y_i(i)) * sig_eps_y +
        eps_tr(trans_i(i)) * sig_eps_tr +
        eps_yt(trans_i(i),y_i(i)) * sig_eps_yt;
        c_i(i) = rgamma(pow(theta,-2), pow(theta,2)*exp(log_chat_i(i)));  // Simulate response
      }
    }
    REPORT(c_i);          // Report the simulation
    REPORT(s_i);          // Stations
    REPORT( s_i_re );       	// Stations with positive counts
    REPORT( y_i );       	// year index
    REPORT(transDist);
    REPORT(trans_i);
    REPORT(eps_y);
    REPORT(eps_s);
    REPORT(eps_tr);
    REPORT(eps_yt);
  }
  
  
  //
  REPORT(log_chat_i);
  REPORT(jnll);
  REPORT(n_i)
  REPORT(p);
  REPORT(eps_s);
  REPORT(sig_eps_s);
  REPORT(eps_y);
  REPORT(sig_eps_y);
  REPORT(rho_y);
  REPORT(eps_tr);
  REPORT(sig_eps_tr);
  REPORT(rho_tr);
  REPORT(eps_yt);
  REPORT(sig_eps_yt);
  REPORT(rho_yt_1);
  REPORT(rho_yt_2);
  REPORT(theta);
  REPORT(ln_lam);
  REPORT(Sigma);  
  REPORT(a);
  return jnll;
}
