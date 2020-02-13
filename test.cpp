// Shows use of vector, matrix and array operations.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data objects
  DATA_INTEGER(i);		// Scalar integer
  
  // Parameter objects
  PARAMETER(p)
  
  // Objects of double type so that they can be printed
  vector<double> v1(1); 
  matrix<double> m1(2,2); 
  m1 << 1,2,3,4;
  vector<double> a1(2); 
  a1 << 11,12;

  v1 = m1 * a1;
  std::cout<<v1<<"\n";
  
  REPORT(v1);
  
  Type ans;
  return ans;
  
  }
