
// create an array with linear spacing
var linspace = function(a, b, n) {
  if (typeof n === "undefined") n = Math.max(Math.round(b - a) + 1, 1);
  if (n < 2) {
    return n === 1 ? [a] : [];
  }
  var i, ret = Array(n);
  n--;
  for (i = n; i >= 0; i--) {
    ret[i] = Math.round((i * b + (n - i) * a) / n * 100) / 100;
  }
  return ret;
};

var randn = function (N1, N2) {
    var u = 0, v = 0;
    while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
}

var ones = function (N1, N2) {
  return 1
}

var add_noise = function (s_bar, tau, the, gam, N_shots, n) {

  N_s = 1 // s_bar.length // length of each signal

  s = s_bar

  for (ii=0; ii<s_bar.length; ii++) {
    // Setup for random variables
    n_P = randn(N_s, N_shots); // standard normal random vector, realizes Poisson noise
    n_G = randn(N_s, N_shots); // standard normal random vector, realizes Gaussian noise

    // Generate observed signals by adding error terms
    s[ii] = s_bar[ii] * ones(1,N_shots) +  // expected average signal (s_bar)
      tau * s_bar[ii] * n +  // shot-to-shot error (delta)
      Math.sqrt(the) * Math.sqrt(s_bar[ii]*ones(1,N_shots) + tau*s_bar[ii]*n) * n_P +  // Poisson noise (p)
      gam * n_G; // Gaussian noise (g)
  }

  return s
}
