
% ADD_NOISE  Simulates signal with Poisson, Gaussian, and shot-to-shot errors. 
%  This considers multiple sources of noise: 
% 
%  INPUTS:
%   s_bar - expected mean signal
%   tau - shot-to-shot std. dev.
%   theta - amplification / scaling factor
%   gamma - Gaussian noise level
%   n_shots - number of signals to generate
% 
%  OUTPUTS:
%   s - set of corrupted signals, with error added
%   Ls - matrix square root of inverse covariance
%   out.s_ave - average of observed signals at each time
%       s_std - standard deviation of observed signals at each time
%       s_tilda - set of single-shot signals (with shot-to-shot error)
% 
%  ------------------------------------------------------------------------
%  
%  POISSON NOISE: 
%   Poisson noise is approximated as Gaussian, with a standard
%   deviation of sqrt of the number of counts. If the data is unscaled
%   the = 1. If the data is scaled, which generally improves
%   inference, theta accounts for the normalization factor. As norm_fact
%   increases, the noise level drops, a consequence of their being more
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

function [s, Ls, out] = ...
    add_noise(s_bar, tau, the, gam, N_shots)

if ~exist('N_shots', 'var'); N_shots = []; end
if isempty(N_shots); N_shots = 1; end % by default only generate one signal


N_s = length(s_bar); % length of each signal


% Setup for random variables
rng(1); % control randomness
n = randn(1, N_shots); % standard normal random variable, realizes shot-to-shot error
n_P = randn(N_s, N_shots); % standard normal random vector, realizes Poisson noise
n_G = randn(N_s, N_shots); % standard normal random vector, realizes Gaussian noise


% Generate observed signals by adding error terms
s = s_bar * ones(1,N_shots) + ... % expected average signal (s_bar)
    tau .* s_bar * n + ... % shot-to-shot error (delta)
    sqrt(the) .* sqrt(s_bar*ones(1,N_shots) + tau.*s_bar*n) .* n_P + ... % Poisson noise (p)
    gam .* n_G; % Gaussian noise (g)


% Calculate statistics of signals
Gs = sparse( ...
    tau^2 + (s_bar*s_bar') + ... % shot-to-shot contribution
    the .* diag(s_bar)) + ... % Poisson contribution
    gam^2 .* speye(N_s);
Ls = chol(inv(Gs)); % get matrix square root

out.s_std = std(s,[],2); % standard deviation over the signals at each time
out.s_ave = mean(s,2); % average of the signals at each time
out.s_tilde = s_bar * ones(1, N_shots)+...
    tau .* s_bar * n; % expected single-shot signal (s_tilda = s_bar + delta)

end

