#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr


rmse  = lambda x,y: np.sqrt(np.mean((x-y)**2, axis = -1))
norm_rmse  = lambda x,y: np.sqrt(np.mean((x-y)**2, axis = -1))/np.std(x, axis = -1)
bias  = lambda x,y: y.mean(axis = -1) - x.mean(axis = -1) 
msess = lambda x,y: 1 - np.sum((x-y)**2, axis = -1)/np.sum((0-x)**2, axis = -1)
mae_wd  = lambda x,y: np.min([np.abs(np.mean(x - y, axis = -1)),
                              np.abs(np.mean((x + 360) - y, axis = -1)),
                              np.abs(np.mean(x - (y + 360), axis = -1))])

bias_wd  = lambda x,y: min([y.mean(axis = -1) - x.mean(axis = -1), 
                               y.mean(axis = -1) - (x.mean(axis = -1) + 360), 
                               (y.mean(axis = -1) + 360) - x.mean(axis = -1)], key = abs)

def auto_corr(gdata):
    acorr = np.correlate(gdata - np.mean(gdata), 
                        gdata - np.mean(gdata), 
                        "full")[len(gdata)-1:]/np.var(gdata) /len(gdata)
    return(acorr[0:10])

def covariance_gufunc(x, y):
    return ((x - x.mean(axis = -1, keepdims = True))
            * (y - y.mean(axis = -1, keepdims = True))).mean(axis = -1)

def pearson_correlation_gufunc(x, y):
    return covariance_gufunc(x, y) / (x.std(axis=-1) * y.std(axis=-1))

## define wrapper functions for apply_ufunc
def pearson_correlation(x, y, dim):
    return xr.apply_ufunc(
        pearson_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='parallelized',
        vectorize=True,
        output_dtypes=[float])

def mae_calculation_wd(x, y, dim):
    return xr.apply_ufunc(
        mae_wd, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

def rmse_calculation(x, y, dim):
    return xr.apply_ufunc(
        rmse, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

def norm_rmse_calculation(x, y, dim):
    return xr.apply_ufunc(
        norm_rmse, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

def bias_calculation(x, y, dim):
    return xr.apply_ufunc(
        bias, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

def bias_calculation_wd(x, y, dim):
    return xr.apply_ufunc(
        bias_wd, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

def msess_calculation(x, y, dim):
    return xr.apply_ufunc(
        msess, x, y,
        input_core_dims=[[dim], [dim]],
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])

def auto_correlation(x, dim):
    return xr.apply_ufunc(
        auto_corr, x, 
        input_core_dims=[[dim]],
        output_core_dims=[["lags"]],
        output_sizes={"lags" : 10},
        dask='allowed',
        vectorize=True,
        output_dtypes=[float])
