
% GET_NOISEL  Evaluate the noise parameters from a set of signals, fitting in log-space.
%  
%  [TAU,THE,GAM] = get_noisel(S) computes the error model parameter
%  for the signals in S, where the second dimension gives repeated
%  measurements of the signal. 
%  
%  [TAU,THE,GAM] = get_noisel(MU,SIG) take the average signal, MU,
%  and standard deviation of the signals, SIG, and computes the error model
%  parameters.
%  
%  [TAU,THE,GAM,FUN] = get_noisel(MU,SIG) updates 
%  
%  AUTHOR: Timothy Sipkens, 2021-05-27

function [tau, the, gam, fun] = get_noisel(s, sig)

% Use covf function with 'pgm' error model.
if ~exist('sig', 'var')
    [~, xlsq] = covf(s, 'pgm', 2);
else
    [~, xlsq] = covf(s, sig);
end

% Quadratic curve used for fitting.
% See Sipkens et al. for why this works.
fun = @(x, s) x(1) .^ 2 .* (s .^ 2) + ...
    x(2) .* s + x(3) .^ 2;

% Interpret polynomial coefficients.
tau = xlsq(1);
the = xlsq(2);
gam = xlsq(3);

end
