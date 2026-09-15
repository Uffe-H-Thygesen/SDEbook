// Inference in a Cox-Ingersoll-Ross process
//
// dX = - lambda*X*dt + sigmaX*dB
//
// based on discrete observations
//
// Y(i) = X(t(i)) + e(i)
//
// where e(i) is N(0,sigmaY^2)
//
// Latent variables are the states.
//
// We use Euler approximation to evalaute transition densities. The time mesh for this
// discretization is finer than the sample interval, i.e. some (many) states are unobserved.

#include <TMB.hpp>

template<class Type>
Type fstrat( const Type & z, const Type & lambda, const Type & xi, const Type & gamma) { return (( lambda * xi - 0.25*gamma*gamma) * exp(-z)  - lambda  ); }

template<class Type>
Type fito( const Type & z, const Type & lambda, const Type & xi, const Type & gamma) { return (( lambda * xi - 0.5*gamma*gamma) * exp(-z)  - lambda  ); }

template<class Type>
Type fsqrt( const Type & z, const Type & lambda, const Type & xi, const Type & gamma)
{ return ( 0.5 * lambda * xi / z - 0.5* lambda * z - 0.125 * gamma*gamma / z); }

template<class Type>
Type g( const Type & z, const Type & gamma) 
{ return( gamma * exp(-0.5 * z ) ); }

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(tsim);     // Time points where X is simulated 
  DATA_VECTOR(iobs);     // Indeces into tsim where X is observed
  DATA_VECTOR(Y);        // Observations taken. Must have same length as iobs.
  DATA_VECTOR_INDICATOR(keep, Y);  // For one-step predictions
  
  PARAMETER_VECTOR(Z);   // States at tsim. Length = length(tsim)

  PARAMETER(lambda);     
  PARAMETER(xi);         
  PARAMETER(loggamma);      
  PARAMETER(logv);      

  Type gamma=exp(loggamma);
  Type v = exp(logv);
  
  Type ans=0;  // ans will be the resulting likelihood

  vector<Type> dt(tsim.size()-1);

  
  for(int i=0;i<dt.size();i++)   
    dt(i) = tsim(i+1)-tsim(i);

  // Include likelihood contributions from state transitions
  for(int i=0;i<dt.size();i++){
    // Type fm = 0.5 * ( fstrat(Z(i),lambda,xi,gamma) + fstrat(Z(i+1),lambda,xi,gamma));
    // Type gm = 0.5 * ( g(Z(i),gamma) + g(Z(i+1),gamma));

    // ans -= dnorm(Z(i+1),Z(i) + fm*dt(i),gm * sqrt(dt(i)),1);

    Type f = fsqrt(Z(i),lambda,xi,gamma);
    Type gi = 0.5*gamma;
    
    ans -= dnorm(Z(i+1),Z(i) + f*dt(i),gi * sqrt(dt(i)),1);
  }

  // Include likelihood contributions from measurements
  for(int i=0;i<Y.size(); ++i){
    int j=CppAD::Integer(iobs(i));
    ans-= keep(i) * dpois(Y(i),v*Z(j) * Z(j),1);
  }

  return ans;
}
