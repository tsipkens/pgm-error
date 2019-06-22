## MATLAB tools for a general TiRe-LII error model
#### (Wat-TiReLII-general-error-model)

This constitutes a general error model software package distributed in
association with work published in the Applied Optics entitled
*[General error model for analysis of laser-induced incandescence
signals](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-56-30-8436)*
by T. A. Sipkens and coworkers.  

----------------------------------------------------------------------

### Description and instructions

This software is an implementation of the general error model
described in the associated work, which simulates signals according
to the described error model (sufficient for generating Figure 1) and
demonstrates a fitting procedure to infer the error model parameters
from the produced signals.

This software package was developed for use with MATLAB 2016a running
on Windows. Users can execute the main script immediately, provided
all attached files are located in the same directory and that
directory is included in the MATLAB path. One must then enter
`main_simulate_C` in the MATLAB command line.

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

#### Size

9.98 KB

----------------------------------------------------------------------

#### Components

This software distribution includes the following files:

README.md		This file.

main_simulate_C.m 	Primary script which loads signals, calls
			function to simulate signals, and performs
			fitting procedure on the data.

simulate_noise.m 	Takes error model parameters and generates
			multiple realizations of signals

J.mat 			A sample expected mean TiRe-LII signal
			generated for C-N2 using the Michelsen
			model in (Michelsen et al., Appl. Phys. B,
			2007)

----------------------------------------------------------------------

#### Contact information:

The primary author of the code is Timothy A. Sipkens, who can be
emailed at [tsipkens@uwaterloo.ca](mailto:tsipkens@uwaterloo.ca).
Alternatively, one can contact Kyle J. Daun at
[kjdaun@uwaterloo.ca](mailto:kjdaun@uwaterloo.ca) who supported code
development.
