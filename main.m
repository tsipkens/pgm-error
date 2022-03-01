
% Main script used to load simualted incandescence, evaluate 
%  the error model, fit error model parameters, and plot results. 

clear;
close all;
% clc;


% Load simulated carbon incandescence trace.
%  Contains time (t) and incandescence (J) produced by 
%  evaluating the Michelsen model in (Michelsen et al., 
%  Appl. Phys. B, 2007) at a wavelength of 500 nm. 
data = csvread('data/lii.csv', 1, 0);
% data = csvread('data/gaus.csv', 1, 0);  % alt. using Gaussian
t = data(:, 1); % time
J = data(:, 2); % incandescence


% Define error model parameters
tau = 0.2; % shot-to-shot variation as a dimensionless std. dev.
the = 1; % amplification / scaling factor
gam = sqrt(2); % Gaussian noise level, in percent of max, i.e. 15 = 15%


% Generate a set of signals with error
N_shots = 500; % number of shots to simulate
s_bar = J .* the; % expected mean signal
[s, ~, G, out] = add_noise(s_bar, tau, the, gam, N_shots);
    % generate observed signals, with error
    

% Fit error model parameters (and display output)
[tau_e, the_e, gam_e, x_var] = get_noise(s); % fit quadratic to variance
[tau_l, the_l, gam_l] = get_noisel(s);
disp('Error model parameters: '); % display results
disp(' ');
fprintf('         <strong> tau      theta    gamma   </strong>\n')
fprintf('          -----    -----    -----   \n')
fprintf('<strong>True</strong>      %4.3f    %4.3f    %4.3f \n', tau, the, gam);
fprintf('<strong>LSQ(log)</strong>  %4.3f    %4.3f    %4.3f \n', tau_l, the_l, gam_l);
fprintf('<strong>LSQ</strong>       %4.3f    %4.3f    %4.3f \n', tau_e, the_e, gam_e);
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
plot_muvar(s, 0, tau_e, the_e, gam_e);

% plot original model parameters
hold on;
max_plot = xlim;
max_plot = max_plot(2);
fplot(@(x) ...
    gam ^ 2 + ...  % Gaussian
    the .* x + ...  % Poisson
    (tau ^ 2) .* (x.^2), ...  % multiplicative
    '-k', [0,max_plot]);
hold off;
h = gca;
h.Legend.String{end} = 'Truth';



%== FIG 3: Plot of average signal versus variance on log-scale ===========%
figure(3); % plot average of observed signals verses variance and fits
plot_muvar(s, 1, tau_e, the_e, gam_e);

% plot original model parameters
hold on;
max_plot = xlim;
max_plot = max_plot(2);
fplot(@(x) gam^2 + the.*x + (tau^2).*(x.^2), ...
    '-k', [0,max_plot]);
hold off;
h = gca;
h.Legend.String{end} = 'Truth';



%== FIG 4: Plot resultant covariance matrix ==============================%
figure(4);
imagesc(G);


