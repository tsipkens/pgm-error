## Tools for PGS error modeling

**Examining *P*oisson, *G*aussian, and *S*hot-to-shot errors.**

This is a software package, originally distributed in association with [Sipkens et al. (2017)][1], which evaluates a general error model for optical signals, including Poisson-Gaussian noise and changes in the measurement conditions between repeated observations (e.g., between laser shots). Particular focus is placed on time-resolved laser-induced (TiRe-LII) and the shot-to-shot variations in quantities like the laser energy and particle volume fraction in the probe volume. This code replaces an archived version available on [figshare](https://figshare.com/articles/MATLAB_tools_for_a_general_TiRe-LII_error_model/5457253/2). 

Two version of this code are available: 

1. A Matlab version code is available in the upper repository. The Matlab code was developed for use with MATLAB 2016a running on Windows and is demonstrated below. Much of the remainder of this README pertains to this version of the code. 

2. A Javascript version of the code is also available in the `js/` directory. The Javascript also uses d3.js to build a web app that demonstrates the error model, available at https://tsipkens.github.io/pgs-error/. 

## Demonstrating the Matlab code

### Simulated signals

The simplest demonstration of this program starts with a theoretical, noiseless signal, considered here with respect to TiRe-LII signals. For simplicity, a sample set of data was included with this distribution for reference. It can be loaded using:

```Matlab
data = csvread('data/lii.csv');
t = data(:, 1); % time
J = data(:, 2); % incandescence
```

This will load three variables into the workspace: (1) `lambda` contains the wavelength for the loaded signal; (2) `t` contains a vector of sample times for the signal; and, most importantly, (3) `J` contains an incandescence trace, evaluated using the Michelsen model from [Michelsen et al. (2007)][2]. Next, define the relevant error model parameters:

```Matlab
% Define error model parameters
tau = 0.2; % shot-to-shot variation as a dimensionless std. dev.
the = 1; % amplification / scaling factor
gam = sqrt(2); % Gaussian noise level
                 % in percent of max, i.e. 15 = 15%
```

The parameters are generally described in [Sipkens et al. (2017)][1]. Briefly, `tau` controls the amount of variation between observations (e.g., between shots), `the` controls the level of Poisson noise (indirectly, by accounting for scaling of the signal), and `gam` controls the cosntant Gaussian component of the background. Now, simulate 500 laser shots:

```Matlab
N_shots = 500; % number of shots to simulate
[s_bar, ~, out] = J .* the; % expected mean signal
[s, ~, out] = add_noise(s_bar, tau, the, gam, N_shots);
    % generate observed signals, with error
```

The second line creates a signal by scaling by the incandescence by the `the` factor above (necessary as the signal and incandescence are given in different units, requiring a conversion). The third line here is the primary call in this program, taking the generated, noiseless signal and creating 500 noisy signals (including shot-to-shot variations as specified by `tau`).

Finally, let's visualize the results. First, plot the fifth realization of the noisy signal:

```Matlab
figure(1);
plot(t, s(:,5), 'b.'); % plot the first noisy signal
```

Now, add the original, average, noiseless signal:

```Matlab
hold on;
plot(t, s_bar, 'k'); % plot overall mean
hold off;
```

Finally, add the noiseless version of the signal, after accounting for the shot-to-shot variations (e.g. that the volume fraction is much less than average):

```Matlab
hold on;
plot(t, out.s_tilde(:,5), 'b'); % plot noiseless signal after shot-to-shot
hold off;
```

### Experimental signals

One can also analyze experimental signals. We will here demonstrate on noisy signals generated from the above procedure, but a simple substitution for experimental signals has been demonstrated by [Sipkens et al. (2017)][1]. Analysis proceeds simply by fitting a polynomial to the mean and variance of the data. That is, compute the mean and standard deviation of the signals, `s`, generated above:

```Matlab
s_ave = mean(s, 2);
s_std = std(s, [], 2);
```

Alternatively, these quantities are output directly from the `add_noise()` function, as shown above. Now, fit a quadratic polynomial to the data using `get_noise()`:

```Matlab
[tau_e, the_e, gam_e, x_var] = get_noise(s);
```

The degree to which the data prescribes to his simple quadratic structure can be demonstrated by plotting the quadratic fit and the data:

```Matlab
figure(2); % plot average of observed signals verses variance and fits
plot(s_ave, s_std.^2, '.'); % plot observed average verses variance
hold on;
max_plot = the * max(J); % maximum of x-axis in plots

% plot quadratic error model fit to variance
fplot(@(x) tau_e^2 + the_e.*x + gam_e.*(x.^2), ...
    '--k', [0,max_plot]);

hold off;
```

[Sipkens et al. (2017)][1] demonstrated that a variety of TiRe-LII signals will show this quadratic relationship.

### A sample script: main.m

Users can execute the main script immediately, provided all attached files are located in the same directory and that directory is included in the MATLAB path. One must then enter `main` in the MATLAB command line. This code will generate Figure 1 from [Sipkens et al. (2017)][1].

The code will display the source error model parameters chosen to generate the data, as well as parameters fit to the signals generated using the given error model. The code will also output two figures: (i) a plot showing a sample signal, expected mean signal, and single-shot expected signal with time (corresponding to Figure 1 in the associated paper) and (ii) a plot show the average verses the variance of the associated signals and the corresponding error model fit to the data.

Error model parameters can be modified by editing the assignment of tau, theta, and gamma in the `main.m` script.

## Data files

Included data files correspond to:

*lii.csv* - Data for a sample expected mean TiRe-LII signal generated for C-N<sub>2</sub> using the Michelsen model from [Michelsen et al. (2007)][mich].

*gaus.csv* - Data corresponding to a Gaussian distribution (used by default in the web app).

----------------------------------------------------------------------

### Contact information

The primary author of the code is Timothy A. Sipkens, who can be
emailed at [tsipkens@uwaterloo.ca](mailto:tsipkens@uwaterloo.ca).
Alternatively, one can contact Kyle J. Daun at
[kjdaun@uwaterloo.ca](mailto:kjdaun@uwaterloo.ca) who supported code
development.

### How to cite

Users of this work should cite,

> Sipkens, T. A., Hadwin, P. J., Grauer, S. J., & Daun, K. J. (2017). General error model for analysis of laser-induced incandescence signals. Applied optics, 56(30), 8436-8445,

and can consider citing this repository, if this code is used directly.

#### References

[Sipkens, T. A., Hadwin, P. J., Grauer, S. J., & Daun, K. J. (2017). General error model for analysis of laser-induced incandescence signals. Applied optics, 56(30), 8436-8445.][1]

[Michelsen, H. A., Liu, F., Kock, B. F., Bladh, H., Boïarciuc, A., Charwath, M., Dreier, T., Hadef, R., Hofmann, M., Reimann, J., & Will, S. (2007). Modeling laser-induced incandescence of soot: a summary and comparison of LII models. Applied physics B, 87(3), 503-521.][2]


[1]: https://www.osapublishing.org/ao/abstract.cfm?uri=ao-56-30-8436

[2]: https://link.springer.com/article/10.1007/s00340-007-2619-5
