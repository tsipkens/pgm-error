
% PLOT_MUVAR  Plot the signal mean v. variance and model parameters, if relevant
%  
%  AUTHOR: Timothy Sipkens, 2021-05-27


function [] = plot_muvar(s, f_log, tau, the, gam)

%-- Parse inputs ----------%
if ~exist('f_log', 'var'); f_log = []; end
if isempty(f_log); f_log = 0; end
%--------------------------%


% Plot observed average verses variance.
plot(mean(s,2), std(s,[],2).^2, '.');

max_plot = the * max(mean(s,2)); % maximum of x-axis in plots


% Add model to plot.
if exist('gam', 'var')
    hold on;

    % Plot quadratic error model fit to variance.
    fplot(@(x) gam ^ 2 + the .* x + (tau ^ 2) .* (x .^ 2), ...
        '--k', [0, max_plot]);

    % Plot only Poisson-Gaussian component of the error model.
    fplot(@(x) gam ^ 2 + the .* x, ...
        '--k', [0, max_plot]);
    
    hold off;
end


% Annotate plot.
xlim([0,max_plot]);
xlabel('<s> [a.u.]');
ylabel('var(s) [a.u.]');
legend('Observed', 'Fit error model', ...
    'Fit Poisson-Gaussian contribution', ...
    'location', 'northwest');


% Change to log-scale, if flagged.
if f_log
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
end

end
