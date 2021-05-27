
% GET_NOISEL  Evaluate the noise parameters from a set of signals, fitting in log-space.
%  
%  [TAU,THE,GAM] = tools.get_noisel(S) computes the error model parameter
%  for the signals in S, where the second dimension gives repeated
%  measurements of the signal. 
%  
%  [TAU,THE,GAM] = tools.get_noisel(MU,SIG) take the average signal, MU,
%  and standard deviation of the signals, SIG, and computes the error model
%  parameters.
%  
%  [TAU,THE,GAM,FUN] = tools.get_noisel(...) adds an output for the 
%  error model curve, FUN.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2021-05-27

function [tau, the, gam, fun] = get_noisel(s, sig)

% Get average and standard deviation is a single input.
if nargin==1
    sig = std(s, [], 2);  % standard deviation
    s = mean(s, 2);  % average
end
s(isinf(s)) = NaN;  % replace infinite values
sig(isinf(sig)) = NaN;
s(s==0) = NaN;
sig(sig==0) = NaN;


% Fit a quadratic curve.
% See Sipkens et al. for why this works.
fun = @(x, s) x(1) .^ 2 .* (s .^ 2) + ...
    x(2) .* s + x(3) .^ 2;

x0 = [0.1, 1, min(sig(sig~=0))];
opts.Display = 'none';
x_lsq = fmincon( ...
    @(x) nansum((log(fun(x, s)) - 2 .* log(sig)).^2), ...
    x0, [], [], [], [], [0,0,0], [], [], opts);

% Interpret polynomial coefficients.
tau = x_lsq(1);
the = x_lsq(2);
gam = x_lsq(3);

end
