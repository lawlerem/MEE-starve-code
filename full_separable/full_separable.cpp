#include <TMB.hpp>
using namespace density;

template<class Type>
matrix<Type> s_covfun(matrix<Type> ds, Type s_range) {
  matrix<Type> ans = ds;
  for(int i=0; i<ds.rows(); i++) {
    for(int j=0; j<ds.cols(); j++) {
      ans(i,j) = s_range * exp(-ds(i,j)/s_range);
    }
  }
  return ans;
};

template<class Type>
matrix<Type> t_covfun(matrix<Type> dt, Type t_ar1) {
  matrix<Type> ans = dt;
  for(int i=0; i<dt.rows(); i++) {
    for(int j=0; j<dt.cols(); j++) {
      ans(i,j) = pow(t_ar1,dt(i,j));
    }
  }
  return ans;
};

template<class Type>
matrix<Type> full_covfun(matrix<Type> ds, matrix<Type> dt, Type st_sd, Type s_range, Type t_ar1) {
  matrix<Type> ans = pow(st_sd,2)*s_covfun(ds,s_range).array()*t_covfun(dt,t_ar1).array();
  return ans;
}


template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y);
  DATA_MATRIX(ds);
  DATA_MATRIX(dt);

  DATA_MATRIX(pred_ds); // Row i gives the distances from pred_re(i) to re
  DATA_MATRIX(pred_dt); // Row i gives the time difference from pred_re(i) to re

  PARAMETER(mu);
  PARAMETER(working_st_sd);
  PARAMETER(working_s_range);
  PARAMETER(working_t_ar1);
  PARAMETER_VECTOR(re);
  PARAMETER_VECTOR(pred_re);

  Type st_sd = exp(working_st_sd);
  ADREPORT(st_sd);

  Type s_range = exp(working_s_range);
  ADREPORT(s_range);

  Type t_ar1 = 2*invlogit(working_t_ar1)-1;
  ADREPORT(t_ar1);

  matrix<Type> full_cov = full_covfun(ds,dt,st_sd,s_range,t_ar1);
  matrix<Type> full_prec = atomic::matinv(full_cov);
  MVNORM_t<Type> re_mvn(full_cov);

  Type nll = 0.0;
  nll += re_mvn(re-mu);
  for(int i=0; i<y.size(); i++) {
    nll -= dpois(y(i),exp(re(i)),true);
  }


  matrix<Type> pred_cross_cov = full_covfun(pred_ds,pred_dt,st_sd,s_range,t_ar1);
  
  for(int i=0; i<pred_re.size(); i++) {
    matrix<Type> c_sigma_inv = matrix<Type>(pred_cross_cov.row(i)) * full_prec;
    Type pred_mu = mu + (c_sigma_inv*(re-mu).matrix())(0,0);
    Type pred_var = pow(st_sd,2) - (c_sigma_inv*matrix<Type>(pred_cross_cov.row(i)).transpose())(0,0);
    nll -= dnorm(pred_re(i), pred_mu, sqrt(pred_var), true);
  }

  vector<Type> response = exp(pred_re);
  ADREPORT(response);

  return nll;
}
