
% ADD_NOISE  Simulates signal with Poisson, Gaussian, and multiplicative errors. 
%  This considers multiple sources of noise: 
%  
%  INPUTS:
%   s_bar     Expected mean signal
%   tau       Multiplicative (shot-to-shot) variance
%   theta     Amplification / scaling factor
%   gamma     Gaussian noise level
%   n_shots   Number of signals to generate
%   seed      Seed for random number generator
%   f_pois    Flag whether to sample Poisson noise from a Poisson distribution.
% 
%  OUTPUTS:
%   s         Set of corrupted signals, with error added
%   Ls        Matrix square root of inverse covariance
%   Gs        Data covariance matrix assocaited with an output
%   out.s_ave     Average of observed signals at each time
%       s_std	  Standard deviation of observed signals at each time
%       s_tilda	  Set of single-shot signals (with shot-to-shot error)
%  
%  ------------------------------------------------------------------------
%  
%  POISSON NOISE: 
%   Poisson noise is approximated as Gaussian, with a standard
%   deviation of sqrt of the number of counts. If the data is unscaled
%   the = 1. If the data is scaled, which generally improves
%   inference, theta accounts for the normalization factor. As norm_fact
%   increases, the noise level drops, a consequence of there being more
%   counts in the data than the scaled b0 indicates.
% 
%  ADDITIVE GAUSSIAN NOISE: 
%   Take additive Gaussian noise as (gam*100)% of peak signal.
%   This represents the minimum noise level and is to model
%   background source, such as electronic noise in the CPC
%   and fluctutions in the aerosol.
%   As a second note, increasing this quantity tends
%   to result in poor reconstruction in background but better
%   reconstruction of the peak of the distribution.
% 
%  AUTHOR: Timothy Sipkens

function [s, L, G, out] = ...
    add_noise(s_bar, tau, the, gam, N_shots, seed, f_pois)

%-- Parse inputs ----------------------------------------%
if ~exist('N_shots', 'var'); N_shots = []; end
if isempty(N_shots); N_shots = 1; end  % by default only generate one signal

if ~exist('seed', 'var'); seed = []; end
if isempty(seed); seed = 1; end

if ~exist('f_pois', 'var'); f_pois = []; end
if isempty(f_pois); f_pois = 0; end  % by defaul, approximate as Gaussian
%--------------------------------------------------------%


N_s = length(s_bar);  % length of each signal


% Setup for random variables
rng(seed);  % control randomness
n = randn(1, N_shots);  % standard normal random variable, realizes shot-to-shot error
n_G = randn(N_s, N_shots);  % standard normal random vector, realizes Gaussian noise
if f_pois == 1  % sample from Poisson distribution
    n_P = (poissrnd(s_bar * ones(1, N_shots)) - s_bar) ./ sqrt(s_bar);  % standardized Poisson r.v.
else  % sample Poisson noise as Gaussian
    n_P = randn(N_s, N_shots);  % standard normal random vector, realizes Poisson noise
end


% Check and replace negative signals.
% Comment this block if negative signals are allowed.
% This may cause a mismatch between Gs and the sampled signals.
f_neg = (tau .* n <= -1);
while any(f_neg)
    n(f_neg) = randn(1, sum(f_neg));
    f_neg = (tau .* n <= -1);
end

% Generate observed signals by adding error terms
s = s_bar .* (1 + ...  % expected average signal (s_bar)
    tau .* n) + ...  % shot-to-shot error (delta)
    sqrt(the) .* sqrt(s_bar .* ones(1,N_shots) + tau .* s_bar .* n) .* n_P + ...  % Poisson noise (p)
    gam .* n_G;  % Gaussian noise (g)


% Calculate statistics of signals
[G, L] = param2cov(tau, the, gam, s_bar);

% Realized, standard deviation, average, and single shot.
out.s_std = std(s,[],2);  % standard deviation over the signals at each time
out.s_ave = mean(s,2);  % average of the signals at each time
out.s_tilde = s_bar * ones(1, N_shots)+...
    tau .* s_bar * n;  % expected single-shot signal (s_tilda = s_bar + delta)

% Expected variance (not realized).
out.s_var = tau ^ 2 .* s_bar .^ 2 + ...
    the .* s_bar + gam ^ 2;


end
