
% PARAM2COV  Get covariance matrix for a set of PGM error model parameters.
%  Build the covariance matrix based on the provided error model
%  parameters. 
%  
%  AUTHOR: Timothy Sipkens, 2021-06-30

function [G, L] = param2cov(tau, the, gam, s_bar)

G = sparse( ...
    tau^2 .* (s_bar * s_bar') + ...  % shot-to-shot contribution
    the .* diag(s_bar)) + ...  % Poisson contribution
    gam^2 .* speye(length(s_bar));

% Get Cholesky factorization of inverse covariance matrix.
% This can be used to weight data in a weighted least-squares scenerio.
if nargout>=2
    L = chol(inv(G));  % get matrix square root
end

end

