
% COV2PARAM  Get covariance matrix for a set of PGM error model parameters.
%  Build the covariance matrix based on the provided error model
%  parameters. 
%  
%  AUTHOR: Timothy Sipkens, 2021-06-30

function [tau, the, gam] = cov2param(G, s)

sig = sqrt(diag(G));  % extract standard deviation

[tau, the, gam] = get_noisel(s, sig);  % get noise parameters

end

