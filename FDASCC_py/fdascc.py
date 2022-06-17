import warnings
warnings.simplefilter('ignore')

from FDASCC_py.setup import setup_fdascc #set_r_home 
setup_fdascc()

import json
import math
import pandas as pd
import numpy as np
from pylab import rcParams
import matplotlib.pyplot as plt

import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import rpy2.robjects.packages as rpackages
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.functions import SignatureTranslatedFunction as STM

from IPython.display import display, Markdown, Latex

def check_fdascc_install():
    """
    Check if FDASCC is installed
    """
    utils = importr('utils')
    installed_packages = utils.installed_packages()
    print(f"FDASCC Installed: {rpackages.isinstalled('FDASCC')}")
    
## available functions
def list_available_fdascc_methods():
    """
    List the methods available in the FDASCC package 
    """
    fdascc = importr('FDASCC')
    methods = list(fdascc.__dict__.keys())[15:24]
    print(methods)
    
# 1d
def scc_1D_one_sample(Y, X, X_band, generate_plots=True):
    try:
        global_env = ro.r.globalenv()
        global_env.clear()
        
        fdascc = importr('FDASCC')
        output = fdascc.scc_1D(Ya=Y, X=X, X_band=X_band)
        # print(Y.shape)
    except:
        return "FDASCC Could not be imported - check install!"
    
    # output data
    
    # out = global_env['out']
    Yhat = output.rx2['Yhat']
    Yhat_pred = output.rx2['Yhat.pred']
    Yhat_deriv = output.rx2['Yhat.deriv']
    Yhat_deriv_pred = output.rx2['Yhat.deriv.pred']
    scc = output.rx2['scc']
    scc_deriv = output.rx2['scc.deriv']
    sce = output.rx2['sce']
    bw = output.rx2['bw']
    bw_deriv = output.rx2['bw.deriv']
    Ya = output.rx2['Ya']
    d_est = output.rx2['d.est']
    d_cov = output.rx2['d.cov']
    derivs = output.rx2['derivs']
    alpha = output.rx2['alpha']
    
    #generate plots
    if generate_plots:
        rcParams['figure.figsize'] = 12,5
        # plot 1: Estimated functions/mean function
        for i in range(1, len(Yhat[0])):
            Yhat_temp = Yhat[0]
            plt.plot(X, Yhat_temp)
            plt.ylabel('fitted function')
            plt.xlabel('t')
            mhat_pred = Yhat_pred
            nalpha = scc.shape[2]
        plt.show()

        #plot 2: scc
        for i in range(nalpha):
            scc_temp = scc[:, :, i]
            plt.plot(X_band, mhat_pred[0])
            plt.ylim(min(scc_temp[:, 0])-1, max(scc_temp[:, 1])+1)
            plt.title(r'Lower SCC when $\alpha$ = {0}'.format(alpha[i]))

            plt.plot(X_band, scc_temp[:, 0], color='blue', linewidth=0.5)
            plt.plot(X_band, scc_temp[:, 1], color='blue', linewidth=0.5)
            plt.xlabel('SCC')
            plt.ylabel('X_Band')
            plt.show()

        #plot 3: derivs
        Yhat_deriv_temp = Yhat_deriv[0]
        plt.plot(X, Yhat_deriv_temp)
        plt.ylabel('fitted derivative')
        plt.xlabel('t')
        plt.show()

        #plot 4: scc for functional derivatives
        if not type(scc_deriv) == rpy2.rinterface_lib.sexp.NULLType:
            mhat_deriv_pred_temp = Yhat_deriv_pred
            for i in range(nalpha):
                scc_deriv_temp = scc_deriv[:, :, i]
                plt.plot(X_band, mhat_deriv_pred_temp[0])
                plt.ylim(min(scc_deriv_temp[:, 0])-1, max(scc_deriv_temp[:, 1])+1)

                plt.plot(X_band, scc_deriv_temp[:,0], color='blue')
                plt.plot(X_band, scc_deriv_temp[:,1], color='blue')
                plt.ylabel("mhat_deriv_pred")
                plt.xlabel("X_band")
                plt.title(r'Lower SCC when $\alpha$ = {0}'.format(alpha[i]))
                plt.show()
        
    
    # print(Yhat.shape)
    return (Yhat, Yhat_pred, Yhat_deriv, Yhat_deriv_pred, scc, scc_deriv, sce, bw, bw_deriv, Ya, d_est, d_cov, derivs, alpha)

# 1d
def scc_1D_two_sample(Ya, Yb, X, X_band, generate_plots=True):
    try:
        global_env = ro.r.globalenv()
        global_env.clear()
        
        fdascc = importr('FDASCC')
        output = fdascc.scc_1D(Ya, Yb, X=X, X_band=X_band)
        # print(Y.shape)
    except:
        return "FDASCC Could not be imported - check install!"
    
    # output data
    
    # out = global_env['out']
    Yhat = output.rx2['Yhat']
    Yhat_pred = output.rx2['Yhat.pred']
    Yhat_deriv = output.rx2['Yhat.deriv']
    Yhat_deriv_pred = output.rx2['Yhat.deriv.pred']
    scc = output.rx2['scc']
    scc_deriv = output.rx2['scc.deriv']
    sce = output.rx2['sce']
    cover_zero = output.rx2['cover.zero']
    bw = output.rx2['bw']
    bw_deriv = output.rx2['bw.deriv']
    knots_est_a = output.rx2['knots.est.a']
    knots_est_b = output.rx2['knots.est.b']
    knots_cov_a = output.rx2['knots.cov.a']
    knots_cov_b = output.rx2['knots.cov.b']
    Ya = output.rx2['Ya']
    Yb = output.rx2['Yb']
    d_est = output.rx2['d.est']
    d_cov = output.rx2['d.cov']
    derivs = output.rx2['derivs']
    alpha = output.rx2['alpha']
    
    #generate plots
    if generate_plots:
        rcParams['figure.figsize'] = 12,5
        # plot 1: Estimated functions/mean function
        for i in range(1, len(Yhat[0])):
            Yhat_temp = Yhat[0]
            plt.plot(X, Yhat_temp)
            plt.ylabel('fitted function')
            plt.xlabel('t')
            mhat_pred = Yhat_pred
            nalpha = scc.shape[2]
        plt.show()

        #plot 2: scc
        for i in range(nalpha):
            scc_temp = scc[:, :, i]
            plt.plot(X_band, mhat_pred[0])
            plt.ylim(min(scc_temp[:, 0])-1, max(scc_temp[:, 1])+1)
            plt.title(r'Lower SCC when $\alpha$ = {0}'.format(alpha[i]))

            plt.plot(X_band, scc_temp[:, 0], color='blue', linewidth=0.5)
            plt.plot(X_band, scc_temp[:, 1], color='blue', linewidth=0.5)
            plt.xlabel('SCC')
            plt.ylabel('X_Band')
            plt.show()

        #plot 3: derivs
        Yhat_deriv_temp = Yhat_deriv[0]
        plt.plot(X, Yhat_deriv_temp)
        plt.ylabel('fitted derivative')
        plt.xlabel('t')
        plt.show()

        #plot 4: scc for functional derivatives
        if not type(scc_deriv) == rpy2.rinterface_lib.sexp.NULLType:
            mhat_deriv_pred_temp = Yhat_deriv_pred
            for i in range(nalpha):
                scc_deriv_temp = scc_deriv[:, :, i]
                plt.plot(X_band, mhat_deriv_pred_temp[0])
                plt.ylim(min(scc_deriv_temp[:, 0])-1, max(scc_deriv_temp[:, 1])+1)

                plt.plot(X_band, scc_deriv_temp[:,0], color='blue')
                plt.plot(X_band, scc_deriv_temp[:,1], color='blue')
                plt.ylabel("mhat_deriv_pred")
                plt.xlabel("X_band")
                plt.title(r'Lower SCC when $\alpha$ = {0}'.format(alpha[i]))
                plt.show()
        
    
    # print(Yhat.shape)
    return (Yhat, Yhat_pred, Yhat_deriv, Yhat_deriv_pred, scc, scc_deriv, sce, cover_zero, bw, bw_deriv, knots_est_a, knots_est_b, 
            knots_cov_a, knots_cov_b, Ya, Yb, d_est, d_cov, derivs, alpha)