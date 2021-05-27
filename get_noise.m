
% GET_NOISE  Evaluate the noise parameters from a set of signals.
%  
%  [TAU,THE,GAM] = tools.get_noise(S) computes the error model parameter
%  for the signals in S, where the second dimension gives repeated
%  measurements of the signal. 
%  
%  [TAU,THE,GAM] = tools.get_noise(MU,SIG) take the average signal, MU,
%  and standard deviation of the signals, SIG, and computes the error model
%  parameters.
%  
%  AUTHOR: Timothy Sipkens

function [tau, the, gam, x_var] = get_noise(s, sig)

% Get average and standard deviation is a single input.
if nargin==1
    sig = std(s, [], 2); % standard deviation
    s = mean(s, 2); % average
end


% Fit a quadratic curve.
% See Sipkens et al. for why this works.
[x_lsq, x_var] = polyfit(s, sig.^2, 2); % fit quadratic to variance


% Interpret polynomial coefficients.
tau = sqrt(x_lsq(1));
the = x_lsq(2);
gam = sqrt(x_lsq(3));


end

