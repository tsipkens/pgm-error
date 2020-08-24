## MATLAB tools for a general optical signal error model

This is a software package, originally distributed in association with [Sipkens et al. (2017)][1], which evaluates a general error model for optical signals, including Poisson-Gaussian noise and changes in the measurement condtitions between repeated observations (e.g. between laser shots). Particular focus is placed on time-resolved laser-induced (TiRe-LII) and the shot-to-shot variations in quantities like the laser energy and particle volume fraction in the probe volume. 

This code replaces an archived version available on [figshare](https://figshare.com/articles/MATLAB_tools_for_a_general_TiRe-LII_error_model/5457253/2). This software package was developed for use with MATLAB 2016a running on Windows. 

### Use

The simplest demonstration of this program starts with a theoretical, noiseless signal, considered here with respect to TiRe-LII signals. For simplicity, a sample set of data was included with this distribution for reference. It can be loaded using:

```Matlab
load('J.mat');
```

This will load three variables into the workspace: (*1*) `lambda` contains the wavelength for the loaded signal; (*2*) `t` contains a vector of sample times for the signal; and, most importantly, (*3*) `J` contains an incandescence trace, evaluated using the Michelsen model from [Michelsen et al. (2007)][2]. Next, define the relevant error model parameters: 

```Matlab
% Define error model parameters
tau = 0.2; % shot-to-shot variation as a dimensionless std. dev.
theta = 1; % amplification / scaling factor
gamma = sqrt(2); % Gaussian noise level
                 % in percent of max, i.e. 15 = 15%
```

The parameters are generally described in [Sipkens et al. (2017)][1]. Briefly, `tau` controls the amount of variation between observations (e.g., between shots), `theta` controls the level of Poisson noise (indirectly, by accounting for scaing of the signal), and `gamma` controls the cosntant Gaussian component of the background. Now, simulate 500 laser shots: 

```Matlab
nn = 500; % number of shots to simulate
s_bar = J.*theta; % expected mean signal
[s,s_ave,s_std,s_tilde] = ...
   simulate_noise(s_bar,tau,theta,gamma,nn);
```

The second line creates a signal by scaling by the incandescence by the `theta` factor above (necessary as the signal and incandescence are given in different units, requiring a conversion). The third line here is the primary call in this program, taking the generated, noiseless signal and creating 500 noisy signals (including shot-to-shot variations as specified by `tau`). 

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
plot(t, s_tilde(:,5), 'b'); % plot noiseless signal after shot-to-shot
hold off;
```

### Components

This software distribution includes the following files:

*README.md* -		This file.

*main_simulate_C.m* - 	Primary script which loads signals, calls
			function to simulate signals, and performs
			fitting procedure on the data.

*simulate_noise.m* -  	Takes error model parameters and generates
			multiple realizations of signals

*J.mat* - 		Matlab data for a sample expected mean TiRe-LII signal
			generated for C-N2 using the Michelsen
			model from [Michelsen et al. (2007)][mich].

### The sample script

Users can execute the main script immediately, provided all attached files are located in the same directory and that directory is included in the MATLAB path. One must then enter `main_simulate_C` in the MATLAB command line. This code will generate Figure 1 from [Sipkens et al. (2017)][1]. 

The code will display the source error model parameters chosen to
generate the data, as well as parameters fit to the signals generated
using the given error model. The code will also output two figures:
(i) a plot showing a sample signal, expected mean signal, and
single-shot expected signal with time (corresponding to Figure 1 in
the associated paper) and (ii) a plot show the average verses the
variance of the associated signals and the corresponding error
model fit to the data.

Error model parameters can be modified by editing the assignment of
tau, theta, and gamma in the `main_simulate_C.m` script.

----------------------------------------------------------------------

#### Contact information

The primary author of the code is Timothy A. Sipkens, who can be
emailed at [tsipkens@uwaterloo.ca](mailto:tsipkens@uwaterloo.ca).
Alternatively, one can contact Kyle J. Daun at
[kjdaun@uwaterloo.ca](mailto:kjdaun@uwaterloo.ca) who supported code
development.

#### How to cite

Users of this work should cite, 

> [Sipkens, T. A., Hadwin, P. J., Grauer, S. J., & Daun, K. J. (2017). General error model for analysis of laser-induced incandescence signals. Applied optics, 56(30), 8436-8445][1], 

and can consider citing this repository, if this code is used directly. 

#### References

[Sipkens, T. A., Hadwin, P. J., Grauer, S. J., & Daun, K. J. (2017). General error model for analysis of laser-induced incandescence signals. Applied optics, 56(30), 8436-8445.][1]

[Michelsen, H. A., Liu, F., Kock, B. F., Bladh, H., Boïarciuc, A., Charwath, M., Dreier, T., Hadef, R., Hofmann, M., Reimann, J., & Will, S. (2007). Modeling laser-induced incandescence of soot: a summary and comparison of LII models. Applied physics B, 87(3), 503-521.][2]


[1]: https://www.osapublishing.org/ao/abstract.cfm?uri=ao-56-30-8436

[2]: https://link.springer.com/article/10.1007/s00340-007-2619-5


```

```