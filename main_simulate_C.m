
% *************************************************************************
% Main script used to load simualted incandescence, evaluate error model, 
%  fit error model parameters, and plot results. 
% *************************************************************************

clear;
close all;
clc;


% Load simulated carbon incandescence trace.
%  Contains incandescence (J), time vector (t), and wavelength (l) produced
%  by evaluating the Michelsen model in (Michelsen et al., 
%  Appl. Phys. B, 2007)
load('J.mat'); 


% Define error model parameters
tau = 0.2; % shot-to-shot variation as a dimensionless std. dev.
theta = 1; % amplification / scaling factor
gamma = sqrt(2); % Gaussian noise level, in percent of max, i.e. 15 = 15%


% Generate a set of signals with error
nn = 500; % number of shots to simulate
s_bar = J.*theta; % expected mean signal
[s,s_ave,s_std,s_tilde] = simulate_noise(s_bar,tau,theta,gamma,nn);
    % generate observed signals, with error
    

% Fit error model parameters
[x_lsq,x_var] = polyfit(s_ave,s_std.^2,2); % fit quadratic to variance
disp('Error model parameters: '); % display results
disp(' ');
fprintf('        tau^2   theta   gamma^2 \n')
fprintf('True    %4.3f   %4.3f    %4.3f \n',tau^2,theta,gamma^3)
fprintf('LSQ     %4.3f   %4.3f    %4.3f \n',x_lsq(1),x_lsq(2),x_lsq(3));
disp(' ');


% Plotting capabilities
figure(1); % plot instances of the signals generated above
[~,ll] = min(s_tilde(1,:)); % signal index to plot, use lowest instance
plot(t,s_bar,'k'); % plots expected mean signal
hold on;
plot(t,s_ave,'.','Color',[0.267,0.005,0.329]);
    % plots average of observed signals
plot(t,s_tilde(:,ll),'--','Color',[0.267,0.6836,0.1328]);
    % plot a set of expected single-shot signals
plot(t,s(:,ll),'.','Color',[0.4531,0.6836,0.1328]);
    % plots a set of realized signals
hold off;
xlabel('t [ns]');
ylabel('s, s_{bar}, s_{tilde}, s_{ave}');
ylim([-20,120]);
legend('s_{bar}','s_{ave}','s_{tilde}','s');


figure(2); % plot average of observed signals verses variance and fits
plot(s_ave,s_std.^2,'.'); % plot observed average verses variance
hold on;
max_plot = theta*max(J); % maximum of x-axis in plots
fplot(@(x) gamma^2+theta.*x+(tau^2).*(x.^2),'-k',[0,max_plot]);
    % plot quadratic error model fit to variance
fplot(@(x) gamma^2+theta.*x,'--k',[0,max_plot]);
    % plot only Poisson-Gaussian component of the error model fit
hold off;
xlim([0,max_plot]);
xlabel('<s> [a.u.]');
ylabel('var(s) [a.u.]');
legend('Observed','General error model','Poisson-Gaussian component',...
    'location','northwest');
