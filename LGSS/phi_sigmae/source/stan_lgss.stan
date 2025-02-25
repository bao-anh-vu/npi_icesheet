data {
  int<lower=0> Tfin;   // # time points (equally spaced)
  vector[Tfin] y;      // log-squared, mean-removed series
  int use_arctanh;
    real sigma_eta; 
    real sigma_eps;
}
parameters {
  real theta_phi;
  vector[Tfin] x;                 // log volatility at time t
}
transformed parameters {
  real<lower = -1, upper = 1> phi;
  
  if (use_arctanh == 1) {
    phi = tanh(theta_phi);
  } else { // use inv logit
    phi = exp(theta_phi) / (1 + exp(-theta_phi));
  }
  
}
model {
  theta_phi ~ normal(0, 1);
  
  x[1] ~ normal(0, sqrt(sigma_eta^2 / (1 - phi^2)));
  for (t in 2:Tfin)
    x[t] ~ normal(phi * x[t - 1], sigma_eta);
  for (t in 1:Tfin)
    y[t] ~ normal(x[t], sigma_eps);
  
}