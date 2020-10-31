
% *************************************************************************
% Main script used to load simualted incandescence, evaluate error model, 
%  fit error model parameters, and plot results. 
% *************************************************************************

clear;
close all;
clc;


% Load simulated carbon incandescence trace.
%  Contains time (t) and incandescence (J) produced by 
%  evaluating the Michelsen model in (Michelsen et al., 
%  Appl. Phys. B, 2007) at a wavelength of 500 nm. 
data = csvread('data_gaus.csv', 1, 0);
t = data(:, 1); % time
J = data(:, 2); % incandescence


% Define error model parameters
tau = 0.2; % shot-to-shot variation as a dimensionless std. dev.
the = 1; % amplification / scaling factor
gam = sqrt(2); % Gaussian noise level, in percent of max, i.e. 15 = 15%


% Generate a set of signals with error
n_shots = 500; % number of shots to simulate
s_bar = J .* the; % expected mean signal
[s, ~, out] = add_noise(s_bar, tau, the, gam, n_shots);
    % generate observed signals, with error
    

% Fit error model parameters (and display output)
[tau_e, the_e, gam_e, x_var] = get_noise(s); % fit quadratic to variance
disp('Error model parameters: '); % display results
disp(' ');
fprintf('        tau     theta   gamma   \n')
fprintf('True    %4.3f   %4.3f    %4.3f \n', tau, the, gam)
fprintf('LSQ     %4.3f   %4.3f    %4.3f \n', tau_e, the_e, gam_e);
disp(' ');



%== FIG 1: A sample LII signal ===========================================%
figure(1); % plot instances of the signals generated above
[~,ll] = min(out.s_tilde(1,:)); % signal index to plot, use lowest instance
plot(t, s_bar, 'k'); % plots expected mean signal
hold on;
plot(t, out.s_ave, '.', 'Color', [0.267,0.005,0.329]);
    % plots average of observed signals
plot(t, out.s_tilde(:,ll), '--', 'Color', [0.267,0.6836,0.1328]);
    % plot a set of expected single-shot signals
plot(t, s(:,ll), '.', 'Color', [0.4531,0.6836,0.1328]);
    % plots a set of realized signals
hold off;
xlabel('t [ns]');
ylabel('s, s_{bar}, s_{tilde}, s_{ave}');
ylim([-20,120]);
legend('s_{bar}', 's_{ave}', 's_{tilde}', 's');



%== FIG 2: Plot of average signal versus variance ========================%
figure(2); % plot average of observed signals verses variance and fits
plot(mean(s,2), std(s,[],2).^2, '.'); % plot observed average verses variance
hold on;
max_plot = the * max(J); % maximum of x-axis in plots

% plot quadratic error model fit to variance
fplot(@(x) gam_e^2 + the_e.*x + (tau_e^2).*(x.^2), ...
    '--k', [0,max_plot]);

% plot original model parameters
fplot(@(x) gam^2 + the.*x + (tau^2).*(x.^2), ...
    '-k', [0,max_plot]);

% plot only Poisson-Gaussian component of the error model
fplot(@(x) gam^2 + the.*x, ...
    '--k', [0,max_plot]);
hold off;

xlim([0,max_plot]);
xlabel('<s> [a.u.]');
ylabel('var(s) [a.u.]');
legend('Observed', 'Fit error model', ...
    'General error model', 'Poisson-Gaussian component', ...
    'location', 'northwest');
