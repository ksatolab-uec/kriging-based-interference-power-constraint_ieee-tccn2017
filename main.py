############################################################################
## Demonstration of Kriging-Based Interference Power Constraint
## Koya Sato, Ph.D.
## 
## Related Article:
## K. Sato and T. Fujii, IEEE Trans. Cogn. Commun. Netw., vol.3, no.1, pp.13-25, March 2017.
## https://ieeexplore.ieee.org/abstract/document/7817747 (Open Access)
#############################################################################

##############################################################################
## NOTE:
## This code simplifies the system model and method from the original article.
## The performances of this demonstration may be different in part
## from results in the article.
## For example:
## 1) Single SU transmitter (TCCN: multiple SUs)
## 2) Binning-based semivariogram modeling (TCCN: residual maximum likelihood)
## (We implemented the original code with C)
##############################################################################

# The MIT License (MIT)
#
# Copyright (c) 2022 Koya SATO.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
from kriging import *
import scipy

'''measurement configuration'''
R_MEASURE   = 100.0           #radius of measurement circle[m]
DCOR        = 20.0            #correlation distance [m]
STDEV_PU    = 8.0             #shadowing standard deviation for PU[dB]
STDEV_SU    = STDEV_PU        #for SU [dB]
PUTX_X      = -1000.0         #x coordinate of primary transmitter
PUTX_Y      = R_MEASURE       #y coordinate of primary transmitter
SUTX_X      = 1000.0          #x coordinate of secondary transmitter
SUTX_Y      = R_MEASURE       #y coordinate of secondary transmitter
PTX         = 30.0            #transmission power [dBm]
ETA         = 3.5             #path loss index
RX_X        = R_MEASURE       #x coordinate of target primary receiver
RX_Y        = R_MEASURE       #y coordinate of target primary receiver

'''parameters for semivariogram modeling and Kriging'''
D_MAX       = 2.0 * R_MEASURE #maximum distance in semivariogram modeling
N_SEMIVAR   = 20  #number for averaging empirical semivariograms

def dosim(sir_d_db, pout, n):
    '''get measurement dataset on a circle'''
    x, y    = gen_location_vector_on_circle(n+1, R_MEASURE)
    x[n]    = RX_X ## evaluation location
    y[n]    = RX_Y

    cov     = gen_varcov_matrix(x, y, DCOR, STDEV_PU)
    z       = gen_multivariate_normal(cov)  #correlated shadowing vector[dB]

    d       = distance(PUTX_X, PUTX_Y, x, y)
    l       = pathloss(d, ETA)              #[dB]
    prx     = PTX - l + z                   #received signal power [dBm]

    '''semivariogram modeling'''
    x_train     = x[:n]
    y_train     = y[:n]
    prx_train   = prx[:n]
    data        = np.vstack([x_train, y_train, prx_train]).T
    d_sv, sv    = gen_emprical_semivar(data, D_MAX, N_SEMIVAR)
    param       = fit_semivar(d_sv, sv)

    '''radio map construction based on Ordinary Kriging'''
    mat     = gen_mat_for_kriging(x_train, y_train, prx_train, param[0], param[1], param[2])
    prx_est, kvar = ordinary_kriging(mat, x_train, y_train, prx_train, RX_X, RX_Y, param[0], param[1], param[2])
    ## kvar: Kriging variance
    err     = prx_est - prx[n] # interpolation error

    '''secondary transmission power design'''
    imax_avg = prx_est - sir_d_db + np.sqrt(2.0 * (kvar + STDEV_SU**2)) * scipy.special.erfinv(2.0 * pout - 1.0)
    d_su  = distance(SUTX_X, SUTX_Y, RX_X, RX_Y)
    l_su  = pathloss(d_su, ETA)
    pmax_est = imax_avg + l_su 

    sir = prx[n] - (pmax_est - l_su + np.random.normal(0.0, STDEV_SU))

    return sir, kvar, err, pmax_est

if __name__ == "__main__":
    loop        = 1000
    pout        = 0.10 # target outage probability
    sir_d_db    = 10.0 # target SIR [dB]
    n           = 50   # number of measurement samples

    cnt_outage = 0
    kstdev_avg = 0.0
    rmse       = 0.0
    psu_avg    = 0.0

    for i in range(loop):
        res         = dosim(sir_d_db, pout, n)
        cnt_outage += res[0] < sir_d_db
        kstdev_avg += np.sqrt(res[1])
        rmse       += res[2]**2
        psu_avg    += 10.0**(0.1 * res[3])
    rmse            = np.sqrt(rmse / loop)

    print("Outage Probability:", cnt_outage / loop)
    print("Average Kriging Standard Deviation:", kstdev_avg / loop, "[dB]")
    print("RMSE:", rmse, "[dB]")
    print("Average Secondary Transmission Power:", 10.0*np.log10(psu_avg / loop), "[dBm]")