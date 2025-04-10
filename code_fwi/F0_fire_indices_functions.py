# -*- coding: utf-8 -*-

# this functions come from https://github.com/steidani/FireDanger/tree/main
# you can also directly use this package

# import modules
import math
import numpy as np
import pandas as pd
import xarray as xr
import datetime as datetime
from numba import jit, vectorize
from typing import Optional, Sequence, Union

# ===============================================
# Canadian Forest Fire Weather Index System (FWI)
# ===============================================

# FWI is initialized with some values for FFMC, DMC and DC components. This means that the first values of the series are not reliable,
# until the index is iterated over several time steps and stabilizes (typically a few days suffice).

# Reference: Wang, Anderson and Suddaby, 2015.
day_lengths = np.array(
  [
      [11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10, 11.2, 11.8],
      [10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2],
      12 * [9],
      [7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1, 8.6, 8.1, 7.8],
      [6.5, 7.5, 9, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8, 7, 6],
  ]
)

drying_factors = np.array(
  [
      [6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8],
      12 * [1.39],
      [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6],
  ]
)

@jit
#@vectorize
def _day_length(lat: Union[int, float], mth: int):  # Union[int, float]
  """Return the average day length for a month within latitudinal bounds."""
  if -30 > lat >= -90:
      dl = day_lengths[0, :]
  elif -15 > lat >= -30:
      dl = day_lengths[1, :]
  elif 15 > lat >= -15:
      return 9
  elif 30 > lat >= 15:
      dl = day_lengths[3, :]
  elif 90 >= lat >= 30:
      dl = day_lengths[4, :]
  elif lat > 90 or lat < -90:
      raise ValueError("Invalid lat specified.")
  else:
      raise ValueError
  return dl[mth - 1]

@jit
def _drying_factor(lat: float, mth: int):  
  """Return the day length factor / drying factor."""
  if -15 > lat >= -90:
      dlf = drying_factors[0, :]
  elif 15 > lat >= -15:
      return 1.39
  elif 90 >= lat >= 15:
      dlf = drying_factors[2, :]
  elif lat > 90 or lat < -90:
      raise ValueError("Invalid lat specified.")
  else:
      raise ValueError
  return dlf[mth - 1]


@vectorize
def ffmc(t, p, w, h, ffmc0):
  """Compute the fine fuel moisture code over one timestep (canadian index).

  Parameters
  ----------
  t: array
    Noon temperature [C].
  p : array
    rainfall amount in open over previous 24 hours, at noon [mm].
  w : array
    Noon wind speed [m/s].
  h : array
    Noon relative humidity [%].
  ffmc0 : array
    Previous value of the fine fuel moisture code.
    ffmc_start = 85.0.
  
  Returns
  -------
  array
    Fine fuel moisture code at the current timestep.

  Example
  -------
  ffmc(17,0,6.944,42,85) == 87.69298009277445

  To Do
  -----
  - if ffmc0 is None, then ffmc0 = ffmc0_start
  - snowcover as parameter (0/1). If snowcover == 1, return ffmc = Nan
  """
  # clip humidity to [0,100]
  #h = np.clip(h, 0, 100) # does not work with @vectorize
  #if ffmc0 is None:
  #    ffmc0 = 85.0

  # convert wind speed from m/s to km/h
  w = w * 3.6

  # Eq. 1
  mo = (147.2 * (101.0 - ffmc0)) / (59.5 + ffmc0) 
  
  if p > 0.5:
      rf = p - 0.5  # Eq.2: Rain reduction to allow for loss in overhead canopy ("effective rainfall")
      if mo > 150.0:
          mo = (mo + 42.5 * rf * np.exp(-100.0 / (251.0 - mo)) * (1.0 - np.exp(-6.93 / rf))
          ) + (0.0015 * (mo - 150.0) ** 2) * np.sqrt(rf)
          # Eq.3b
      elif mo <= 150.0:
          mo = mo + 42.5 * rf * np.exp(-100.0 / (251.0 - mo)) * (1.0 - np.exp(-6.93 / rf))
          # Eq.3a: The real moisture content of pine litter ranges up to about 250 percent, so we cap it at 250
      if mo > 250.0:
          mo = 250.0
  # Eq.4: Equilibrium moisture content from drying phase
  ed = (
      0.942 * (h ** 0.679)
      + (11.0 * np.exp((h - 100.0) / 10.0))
      + 0.18 * (21.1 - t) * (1.0 - 1.0 / np.exp(0.1150 * h))
  ) 

  if mo < ed:
        # Eq. 5 Equilibrium moisture content from wetting phase
      ew = (
          0.618 * (h ** 0.753)
          + (10.0 * np.exp((h - 100.0) / 10.0))
          + 0.18 * (21.1 - t) * (1.0 - 1.0 / np.exp(0.115 * h))
      )  
      if mo < ew:
          #Eq. 7a (ko) Log wetting rate at the normal temperature of 21.1 C 
          kl = 0.424 * (1.0 - ((100.0 - h) / 100.0) ** 1.7) + (
              0.0694 * np.sqrt(w)
          ) * (1.0 - ((100.0 - h) / 100.0) ** 8)
          # Eq. 7b Affect of temperature on wetting rate
          kw = kl * (0.581 * np.exp(0.0365 * t))  
          # Eq. 9
          m = ew - (ew - mo) / 10.0 ** kw  
      elif mo > ew:
          m = mo
  elif mo == ed:
      m = mo
  else: # if mo > ed
        #Eq. 6a (ko) Log drying rate at the normal temperature of 21.1 C
      kl = 0.424 * (1.0 - (h / 100.0) ** 1.7) + (0.0694 * np.sqrt(w)) * (
          1.0 - (h / 100.0) ** 8
      )
      # Eq. 6b Affect of temperature on  drying rate
      kw = kl * (0.581 * np.exp(0.0365 * t))
      # Eq.8
      m = ed + (mo - ed) / 10.0 ** kw

  # Eq. 10 Final ffmc calculation
  ffmc = (59.5 * (250.0 - m)) / (147.2 + m)
  # Constraints: ffmc is scaled between 0 and 101
  # ffmc = min(max(0.0,ffmc),101.0)
  if ffmc > 101.0:
      ffmc = 101.0
  elif ffmc <= 0.0:
      ffmc = 0.0
  
  return ffmc


@vectorize
def ffmc2(t, p, w, h, ffmc0):
  """Compute the fine fuel moisture code over one timestep (canadian index).

  Parameters
  ----------
  t: array
    Noon temperature [C].
  p : array
    rainfall amount in open over previous 24 hours, at noon [mm].
  w : array
    Noon wind speed [m/s].
  h : array
    Noon relative humidity [%].
  ffmc0 : array
    Previous value of the fine fuel moisture code.
    ffmc_start = 85.0.
  
  Returns
  -------
  array
    Fine fuel moisture code at the current timestep.

  Example
  -------
  ffmc(17,0,6.944,42,85) == 87.69298009277445

  To Do
  -----
  - if ffmc0 is None, then ffmc0 = ffmc0_start
  - snowcover as parameter (0/1). If snowcover == 1, return ffmc = Nan
  """
  # clip humidity to [0,100]
  #h = np.clip(h, 0, 100) # does not work with @vectorize
  #if ffmc0 is None:
  #    ffmc0 = 85.0

  # convert wind speed from m/s to km/h
  w = w * 3.6

  # Eq. 1
  mo = (147.27723 * (101.0 - ffmc0)) / (59.5 + ffmc0) 
  
  if p > 0.5:
      rf = p - 0.5  # Eq.2: Rain reduction to allow for loss in overhead canopy ("effective rainfall")
      if mo > 150.0:
          mo = (mo + 42.5 * rf * np.exp(-100.0 / (251.0 - mo)) * (1.0 - np.exp(-6.93 / rf))
          ) + (0.0015 * (mo - 150.0) ** 2) * np.sqrt(rf)
          # Eq.3b
      elif mo <= 150.0:
          mo = mo + 42.5 * rf * np.exp(-100.0 / (251.0 - mo)) * (1.0 - np.exp(-6.93 / rf))
          # Eq.3a: The real moisture content of pine litter ranges up to about 250 percent, so we cap it at 250
      if mo > 250.0:
          mo = 250.0
  # Eq.4: Equilibrium moisture content from drying phase
  ed = (
      0.942 * (h ** 0.679)
      + (11.0 * np.exp((h - 100.0) / 10.0))
      + 0.18 * (21.1 - t) * (1.0 - 1.0 / np.exp(0.1150 * h))
  ) 

  if mo < ed:
        # Eq. 5 Equilibrium moisture content from wetting phase
      ew = (
          0.618 * (h ** 0.753)
          + (10.0 * np.exp((h - 100.0) / 10.0))
          + 0.18 * (21.1 - t) * (1.0 - 1.0 / np.exp(0.115 * h))
      )  
      if mo < ew:
          #Eq. 7a (ko) Log wetting rate at the normal temperature of 21.1 C 
          kl = 0.424 * (1.0 - ((100.0 - h) / 100.0) ** 1.7) + (
              0.0694 * np.sqrt(w)
          ) * (1.0 - ((100.0 - h) / 100.0) ** 8)
          # Eq. 7b Affect of temperature on wetting rate
          kw = kl * (0.581 * np.exp(0.0365 * t))  
          # Eq. 9
          m = ew - (ew - mo) / 10.0 ** kw  
      elif mo > ew:
          m = mo
  elif mo == ed:
      m = mo
  else: # if mo > ed
        #Eq. 6a (ko) Log drying rate at the normal temperature of 21.1 C
      kl = 0.424 * (1.0 - (h / 100.0) ** 1.7) + (0.0694 * np.sqrt(w)) * (
          1.0 - (h / 100.0) ** 8
      )
      # Eq. 6b Affect of temperature on  drying rate
      kw = kl * (0.581 * np.exp(0.0365 * t))
      # Eq.8
      m = ed + (mo - ed) / 10.0 ** kw

  # Eq. 10 Final ffmc calculation
  ffmc = (59.5 * (250.0 - m)) / (147.2 + m)
  # Constraints: ffmc is scaled between 0 and 101
  # ffmc = min(max(0.0,ffmc),101.0)
  if ffmc > 101.0:
      ffmc = 101.0
  elif ffmc <= 0.0:
      ffmc = 0.0
  
  return ffmc




@vectorize
def dmc(t, p, h, mth: int, lat: float, dmc0: float): 
  """Compute the Duff moisture code over one time step (canadian index).

  Parameters
  ----------
  t: array
    Noon temperature [C].
  p : array
    rainfall amount in open over previous 24 hours, at noon [mm].
  h : array
    Noon relative humidity [%].
  mth : integer array
    Month of the year [1-12].
  lat : float
    Latitude in degrees.
  dmc0 : float
    Previous value of the Duff moisture code.

  Returns
  -------
  array
    Duff moisture code at the current timestep

  Example
  ------- 
  dmc(17,0,42,6,45.98,6) == 8.5450511359999997
  """
  # clip humidity to [0,100]
  #h = np.clip(h, 0, 100)
  
  #if dmc0 is None:
  #    dmc0 = 6
 #import pdb; pdb.set_trace()
  if np.isnan(p) or np.isnan(t) or np.isnan(h) :
      dmc = np.nan
  else:
      if t < -1.1:
          rk = 0
      else:
          dl = _day_length(lat, mth)
          # Eqs.16 and 17
          rk = 1.894 * (t + 1.1) * (100.0 - h) * dl * 0.0001  
    
      if p > 1.5:
          ra = p
          # Eq.11 Effective rainfall
          rw = 0.92 * ra - 1.27  
          # Eq.12 from R-package cffdrs
          wmi = 20.0 + 280.0 / np.exp(0.023 * dmc0) 
          if dmc0 <= 33.0:
              # Eq.13a
              b = 100.0 / (0.5 + 0.3 * dmc0)  
          else: # dmc0 > 33.0
              if dmc0 <= 65.0:
                  # Eq.13b
                  b = 14.0 - 1.3 * np.log(dmc0)  
              else:
                    # Eq.13c
                  b = 6.2 * np.log(dmc0) - 17.2 
          # Eq.14 duff moisture content after p
          wmr = wmi + (1000 * rw) / (48.77 + b * rw)  
          # Eq.15
          pr = 43.43 * (5.6348 - np.log(wmr - 20.0))  
      else:  # p <= 1.5
          pr = dmc0
      
      if pr < 0.0:
          pr = 0.0
      # Calculate final dmc
      dmc = pr + rk
      # Constraints: dmc is scaled between max(0, dmc)
      if dmc < 0:
          dmc = 0.0
          
  return dmc

@vectorize
def dc(t, p, mth, lat, dc0):  
  """Compute the drought code over one time step (canadian index).

  Parameters
  ----------
  t: array
    Noon temperature [C].
  p : array
    rainfall amount in open over previous 24 hours, at noon [mm].
  mth : integer array
    Month of the year [1-12].
  lat : float
    Latitude.
  dc0 : float
    Previous value of the drought code.

  Returns
  -------
  array
    Drought code at the current timestep

  Example
  ------- 
  dc(17,0,4,45.98,15) == 19.013999999999999
  """
  if np.isnan(p) or np.isnan(t):
      dc = np.nan
  else:
      fl = _drying_factor(lat, mth) # influence of latitude, from R-package cffdrs
    
      if t < -2.8:
          t = -2.8
      # Eq.22 Potential Evapotranspiration
      pe = (0.36 * (t + 2.8) + fl) / 2  
      if pe < 0.0:
          pe = 0.0
    
      if p > 2.8:
          ra = p
          # Eq.18 Effective rainfall
          rw = 0.83 * ra - 1.27  
          # Eq.19 Moisture equivalent of the previous day's DC
          smi = 800.0 * np.exp(-dc0 / 400.0)  
          # Eqs.20 and 21
          dr = dc0 - 400.0 * np.log(1.0 + ((3.937 * rw) / smi))  
          if dr > 0.0:
              dc = dr + pe
          elif np.isnan(dc0):
              dc = np.NaN
          else:
              dc = pe
      else:  # if precip is less than 2.8 then use yesterday's DC
          dc = dc0 + pe
  return dc

def isi(w, ffmc):
  """Initialize spread index (canadian index).

  Parameters
  ----------
  w : array
    Noon wind speed [m/s].
  ffmc : array
    Fine fuel moisture code.

  Returns
  -------
  array
    Initial spread index.

  Example
  ------- 
  isi(6.944444444444445,87.6929800927744) == 10.853661073655068
  """
  # convert wind speed from m/s to km/h
  w = w * 3.6
  # Eq.1  Moisture content
  mo = 147.2 * (101.0 - ffmc) / (59.5 + ffmc)  
  # Eq.25 Fine Fuel Moisture
  ff = 19.1152 * np.exp(mo * -0.1386) * (1.0 + (mo ** 5.31) / 49300000.0)  
  # Eq.26 Spread Index Equation (with Wind Effect)
  isi = ff * np.exp(0.05039 * w)  
  return isi

def bui(dmc, dc):
  """Build-up index (canadian index).

  Parameters
  ----------
  dmc : array
    Duff moisture code.
  dc : array
    Drought code.

  Returns
  -------
  array
    Build up index.
  
  Example
  ------- 
  bui(8.5450511359999997,19.013999999999999) == 8.4904265358371838
  """
  bui = np.where(
      dmc <= 0.4 * dc,
      # Eq.27a
      (0.8 * dc * dmc) / (dmc + 0.4 * dc),  
      # Eq.27b
      dmc - (1.0 - 0.8 * dc / (dmc + 0.4 * dc)) * (0.92 + (0.0114 * dmc) ** 1.7),
  )  
  return np.clip(bui, 0, None)


def fwi(isi, bui):
  """Fire weather index in S-scale (canadian index).

  Parameters
  ----------
  isi : array
    Initial spread index
  bui : array
    Build up index.

  Returns
  -------
  array
    Fire weather index.
  
  Example
  ------- 
  fwi(10.853661073655068,8.4904265358371838) = 10.096371392382368
  """
  fwi = np.where(
      bui <= 80.0,
      # Eq.28a
      0.1 * isi * (0.626 * bui ** 0.809 + 2.0),  
      # Eq.28b
      0.1 * isi * (1000.0 / (25.0 + 108.64 / np.exp(0.023 * bui))),
  )  
  # Eqs.30a and 30b Constraint if fwi > 1
  fwi[fwi > 1] = np.exp(2.72 * (0.434 * np.log(fwi[fwi > 1])) ** 0.647)  
  return fwi

def daily_severity_rating(fwi):
  """Daily severity rating (canadian index).

  Parameters
  ----------
  fwi : array
    Fire weather index

  Returns
  -------
  array
    Daily severity rating.
  """
  return 0.0272 * fwi ** 1.77

