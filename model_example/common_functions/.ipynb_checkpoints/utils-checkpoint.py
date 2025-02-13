import pandas as pd
import numpy as np
import os

from scipy.stats import pearsonr
from scipy.integrate import simpson

from sklearn.pipeline import Pipeline
from sklearn.feature_selection import RFE, SelectKBest
from sklearn.model_selection import KFold, StratifiedGroupKFold, train_test_split, StratifiedKFold, GridSearchCV, cross_val_score, cross_validate, RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import (fbeta_score, make_scorer, 
precision_score, accuracy_score, confusion_matrix, ConfusionMatrixDisplay, RocCurveDisplay, precision_recall_curve,
mean_squared_error, r2_score)
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle

# scikit-optimize bayesian optimization packages
from skopt import BayesSearchCV
from skopt.space import Real, Categorical, Integer

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly
import plotly.io as pio
from plotly.subplots import make_subplots

pio.renderers.default = "svg"

# From an input file (xvg from PyReweighting-1D), ingests the GaMD pmf from PyReweighting as a Pandas dataframe
# 
# file - string, file name containing 1-D pmf data from PyReweighting-1D.py
def get_pmf(file):
    x,y = np.loadtxt(file,comments=["@", "#"],unpack=True)
    heat = pd.DataFrame(np.asarray([x, y]).T, columns=["X", "Y"])
    heat["pop"] = np.exp(-heat["Y"])
    return heat

# From an input file (xvg from PyReweighting-2D), ingests the GaMD pmf from PyReweighting as a Pandas dataframe
# 
# file - string, file name containing 2-D pmf data from PyReweighting-2D.py
def get_pmf_2d(file):
    x,y,z = np.loadtxt(file,comments=["@", "#"],unpack=True)
    heat = pd.DataFrame(np.asarray([x, y, z]).T, columns=["X", "Y", "Z"])
    heat["pop"] = np.exp(-heat["Z"])
    return heat

# Plots 1-D GaMD pmf from ingested pandas dataframe
#
# xlabel, ylabel - plotly axis titles for x and y axes
# scaler - scaling factor for PMF trace to adjust to match histogram height
def plot_pmf(pmfDf, scaler):
    px.line(pmfDf, "X", pmfDf["pop"] * scaler).show()
    
# Plots 1-D GaMD pmf overlaid onto TI histogram for a particular degree of freedom
#
# pmfDf - 1-D pmf output from PyReweighting-1D.py
# ti_df - lambda production data from thermodynamic integration
# dof_str - string of DOF in ti_df that is to be plotted as a histogram
# scaler - scaling factor for PMF trace to adjust to match histogram height
# bins - amount of bins desired in TI DOF histogram
def plot_pmf_TI_1D_fig(pmfDf, ti_df, dof_str, scaler, bins):
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x = pmfDf["X"],
        y = pmfDf["pop"] * scaler,
        mode="lines",
        line=dict(width=4, color="#006d2c"),
        name="-exp(GaMD PMF)"
    ))
    
    hist_data, bin_edges = np.histogram(ti_df[dof_str], bins=bins, density=False)
    total = max(hist_data)
    factor = scaler/total
        
    fig.add_trace(
        go.Bar(
            x=bin_edges[:-1,],
            y=hist_data * factor,
            width=np.diff(bin_edges),
            name="TI Lambda Prod Histogram",
            # histnorm="probability",
            marker=dict(color="#edf8e9", line=dict(color="black", width=1))
        )
    )
    
    return fig


# Plots 2-D GaMD pmf (only GaMD, no TI data) for 2 particular degrees of freedom, formatted pretty
# 
# pmfDf - 2-D pmf output from PyReweighting-2D.py
# xlabel, ylabel - plotly axis titles for x and y axes
def plot_pmf_2d_fig(pmfDf, xlabel, ylabel):
    pivot = pmfDf.pivot("Y", "X", "Z")
    fig = go.Figure(
        go.Contour(
            z = pivot, 
            x=pivot.columns, 
            y=pivot.index,
            colorscale=plotly.colors.sequential.Sunset_r,
            contours=dict(
                start=1,
                end=10,
                size=1,
                
            ),
            line=dict(color="black", width=1),
            showscale=True,
            colorbar=dict(
                outlinecolor="black",
                outlinewidth=2,
                ypad=30,
                tickmode="linear",     
                title=dict(
                    side="bottom",
                    text="Free energy <br>(kcal/mol)</br>", 
                    font=dict(size=32)),
                tickfont=dict(size=32)
            )
        )
    )
    return fig

# Plots 2-D TI data as a scatterplot overlaid onto GaMD pmf for 2 particular degrees of freedom 
#
# pmfs - 2-D pmf output from PyReweighting-2D.py
# ti_df - lambda production data from thermodynamic integration
# dof1, dof2 - desired TI DOFs to be plotted 
def plot_pmf_with_TI_2d_fig(pmfs, ti_df, dof1, dof2):
    pivot = pmfs.pivot("Y", "X", "Z")
    fig = go.Figure(
        go.Contour(
            z = pivot, 
            x=pivot.columns, 
            y=pivot.index,
            colorscale=plotly.colors.sequential.Sunset_r,
            contours=dict(
                start=1,
                end=10,
                size=1,
                
            ),
            line=dict(color="black", width=1),
            showscale=True,
            colorbar=dict(
                outlinecolor="black",
                outlinewidth=2,
                ypad=30,
                tickmode="linear",     
                title=dict(
                    side="bottom",
                    text="Free energy <br>(kcal/mol)</br>", 
                    font=dict(size=32)),
                tickfont=dict(size=32)
            )
        )
    ).add_trace(go.Scatter(
        x=ti_df[dof1], y=ti_df[dof2], mode="markers", 
        marker=dict(
            size=3, color="white", line=dict(color="black", width=0.5))
    ))
    return fig

# Shuffles the TI data, then does a 5x5-fold cross validation (Repeated Stratified k-fold CV) to 
# produce a training and test R2 as well as a feature importance dataframe containing the feature importances
# in each of the 25 separate trainings of the model.
# 
# pipeline - pre-defined model pipeline (we are using recursive feature elimination via decision tree, then a random forest)
# Xin - independent variables, here geometric degrees of freedom (aka interatomic distances and rotamer angles)
# Yin - dependent variables, here weighted DV/DL values weighted by Gaussian quadrature
# lambdas - lambda values for each frame, used for stratification
def benchmark_model(pipeline, Xin, Yin, lambdas):
    X_shf, Y_shf = shuffle(Xin, Yin, random_state=42)
    benchmark = cross_validate(
        pipeline, X_shf, Y_shf,
        cv=RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42).get_n_splits(X_shf, lambdas), 
        n_jobs=-1, verbose=0, return_train_score=True, return_estimator=True)
    print("Avg. training r2: ")
    print(round(np.mean(benchmark["train_score"]), 4))
    print("Training r2 std dev: ")
    print(round(np.std(benchmark["train_score"])/np.sqrt(len(benchmark["train_score"])), 4))
    print("Avg. test r2: ")
    print(round(np.mean(benchmark["test_score"]), 4))
    print("Testing r2 std dev: ")
    print(round(np.std(benchmark["test_score"])/np.sqrt(len(benchmark["test_score"])), 4))

    fi_set = False
    counter = 0
    
    # creating a dataframe for feature importances for each of the 25 separate model training iterations
    for pipe in benchmark["estimator"]:
        if (not fi_set):
            # For the first iteration, sets up the DF, then populates the df with feature importances
            orig_names = pipe.steps[0][1].get_feature_names_out()
            fi = pd.DataFrame(pipe.steps[2][1].feature_importances_)
            fi.index = [orig_names[i] for i in range(len(orig_names)) if pipe.steps[1][1].support_[i]]
            fi_set = True
            counter += 1
        else:
            fi[counter] = pipe.steps[2][1].feature_importances_
            counter+= 1
            
    fi["Mean"] = fi.mean(axis=1)
    fi["Median"] = fi.median(axis=1)
    return fi
    
    


# takes TI data, takes a 2-D GaMD pmf, and returns True/False for each row of TI data whether the 
# data falls within the favorable region in the pmf or not.

# thr - threshold for "favorable" (GaMD will output pmf values, the most favorable area will have pmf = 0, the least favorable
# will be some PMF max value
# xeps and yeps are the grid spacing of the pmf - we reweighted our rotamer distributions in GaMD by discretizations of 6 degrees.
# thus, xeps and yeps are set to 3.
def rot_in_pmf_cont(df, x_dof, y_dof, pmf, thr):
    in_pmf = []
    
    pivot = pmf.pivot("Y", "X", "Z")
    x = pivot.columns
    y = pivot.index
    X, Y = np.meshgrid(x, y)
    
    for i,j in zip(df[x_dof], df[y_dof]):
        
        nearest_x = x[np.argmin(np.abs(x - i))]
        nearest_y = y[np.argmin(np.abs(y - j))]
                        
        if pivot.loc[nearest_y, nearest_x] <= thr:
            in_pmf.append(True)
        else:
            in_pmf.append(False) 
    
    return in_pmf



    

# cd into the directory you want, then find the Gaussian-quadrature weighted lambdas using this function. 

# We are assuming a 12 point Gaussian quadrature 

# rst = restart, which is default, however certain jobs (such as those with NMR restraints added to 
# previously equilibrated coordinates) are not restarted from previous trajectories. One extra
# DV/DL is added for non-restart jobs which must be sliced.
def weighted_lambdas(runNum, rst=True):
    weights = [0.02359, 0.05347, 0.08004, 0.10158, 0.11675, 0.12457, 0.12457, 0.11675, 0.10158, 0.08004, 0.05347, 0.02359]

    dvdls = []
    run_integrals = []
    weighted_lam_list = []
    for run in range(1, runNum + 1):
        # list of 12 lambda averages
        averages = []

        for y in range(12):
            lambda1 = os.popen(f"cat ./run_{run}/{y}_lambda/lambda_prod.out | grep \'DV/DL\'").read().split('\n')
            
            if rst:
                total_lambda = [float(l[9:].strip()) for l in lambda1[:-7]][::2]
            else:
                # non-restart job, first value must be sliced 
                total_lambda = [float(l[9:].strip()) for l in lambda1[2:-7]][::2]
                
            dvdls.append(total_lambda)
            weighted_lams = [weights[y] * i for i in total_lambda]
            weighted_lam_list.append(weighted_lams)
            averages.append(np.mean(total_lambda))

    
    return (weighted_lam_list, dvdls)


# cd into the directory you want, then find the DV/DLs at each checkpoint step using this function

# rst = restart, which is default, however certain jobs (such as those with NMR restraints added to 
# previously equilibrated coordinates) are not restarted from previous trajectories. One extra
# DV/DL is added for non-restart jobs which must be sliced.
def read_rotamers(runNum, rst=True):
    chi = pd.DataFrame()
    
    weighted_lam_list = weighted_lambdas(runNum, rst=rst)
    
    weight_lam_flat = np.concatenate(weighted_lam_list[0]).ravel()
    chi["weight_dvdl"] = weight_lam_flat
    chi["dvdl"] = np.concatenate(weighted_lam_list[1]).ravel()

    runs = [[i] * 2400 for i in range(1, runNum + 1)]
    lambdas = [[j] * 200 for j in range(1, 13)] * (runNum)

    chi["Run"] = np.concatenate(runs).ravel()
    chi["Lambda"] = np.concatenate(lambdas).ravel()
    
    return chi
    




