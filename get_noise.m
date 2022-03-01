
% GET_NOISE  Evaluate the noise parameters from a set of signals.
%  
%  [TAU,THE,GAM] = get_noise(S) computes the error model parameter
%  for the signals in S, where the second dimension gives repeated
%  measurements of the signal. 
%  
%  [TAU,THE,GAM] = get_noise(MU,SIG) take the average signal, MU,
%  and standard deviation of the signals, SIG, and computes the error 
%  model parameters.
%  
%  [TAU,THE,GAM,X_STD] = get_noise(...) add output for uncertainties
%  on error model parameters. 
%  
%  AUTHOR: Timothy Sipkens

function [tau, the, gam, x_std] = get_noise(s, sig)

% Get average and standard deviation is a single input.
if nargin==1
    sig = std(s, [], 2); % standard deviation
    s = mean(s, 2); % average
end

% Fit a quadratic curve.
% See Sipkens et al. for why this works.
[x_lsq, x_std] = polyfit(s, sig.^2, 2); % fit quadratic to variance

% Interpret polynomial coefficients.
tau = sqrt(x_lsq(1));
the = x_lsq(2);
gam = sqrt(x_lsq(3));

%{
% Uncertainties. 
x_std = (inv(x_std.R)*inv(x_std.R)') * ...
    x_std.normr ^ 2 / x_std.df;
x_std = sqrt(diag(x_std));

x_std(1) = std(sqrt( ...
    normrnd(tau, x_std(1), [1e5,1])));
x_std(3) = std(sqrt( ...
    normrnd(gam, x_std(3), [1e5,1])));
%}

end

