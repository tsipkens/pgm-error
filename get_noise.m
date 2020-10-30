
% GET_NOISE  Evaluate the noise parameters from a set of signals.
%=========================================================================%

function [tau, the, gam, x_var] = get_noise(s)

% Get average and standard deviation.
s_ave = mean(s, 2); % average
s_std = std(s, [], 2); % standard deviation


% Fit a quadratic curve.
% See Sipkens et al. for why this works.
[x_lsq, x_var] = polyfit(s_ave, s_std.^2, 2); % fit quadratic to variance


% Interpret polynomial coefficients.
tau = sqrt(x_lsq(1));
the = x_lsq(2);
gam = sqrt(x_lsq(3));


end

