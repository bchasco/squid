#include <TMB.hpp>
// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR (ef);
  DATA_VECTOR (idx);
  DATA_VECTOR (ef_sd);
  DATA_VECTOR (idx_sd);
  
  PARAMETER_VECTOR (u);
  PARAMETER(a);
  PARAMETER(b);
  Type nll=0.;
  
  int ni = ef.size();
  for(int i=0;i<ni;i++){
    nll -= dnorm(u(i),idx(i),idx_sd(i),true);
    nll -= dnorm(ef(i),a+b*u(i),ef_sd(i),true);
  }
  
  vector<Type> p(100);
  for(int i=0;i<100;i++){
    p(i) = a+b*i*0.2;
  }
  REPORT(u);
  REPORT(a);
  REPORT(b);
  ADREPORT(p);
    
  return nll;
  
}
