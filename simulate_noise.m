function [s,s_ave,s_std,s_tilde] = simulate_noise(s_bar,tau,theta,gamma,nn)
% SIMULATE_NOISE Simulates signal with Poisson, Gaussian, and shot-to-shot 
% errors. 

% INPUTS:   s_bar - expected mean signal
%           tau - shot-to-shot std. dev.
%           theta - amplification / scaling factor
%           gamma - Gaussian noise level
%           nn - number of signals to generate

% OUTPUTS:  s - set of observed signal with error added
%           s_ave - average of observed signals at each time
%           s_std - standard deviation of observed signals at each time
%           s_tilda - set of single-shot signals (with shot-to-shot error)

% Setup for random variables
mm = length(s_bar); % length of each signal
n = randn(1,nn); % standard normal random variable, 
    % realizes shot-to-shot error
nP = randn(mm,nn); % standard normal random vector, realizes Poisson noise
nG = randn(mm,nn); % standard normal random vector, realizes Gaussian noise

% Generate observed signals by adding error terms
s = s_bar*ones(1,nn)+... % expected average signal (s_bar)
    tau.*s_bar*n+... % shot-to-shot error (delta)
    sqrt(theta).*sqrt(s_bar*ones(1,nn)+tau.*s_bar*n).*nP+... % Poisson noise (p)
    gamma.*nG; % Gaussian noise (g)

% Calculate statistics of signals
s_std = std(s,[],2); % standard deviation over the signals at each time
s_ave = mean(s,2); % average of the signals at each time
s_tilde = s_bar*ones(1,nn)+...
    tau.*s_bar*n; % expected single-shot signal (s_tilda = s_bar + delta)

end

