'''
Palaeoclimate constraints tighten future Arctic winter sea ice projections
C. Y. Fu, M. B. Osman, M. A. Aquino-López
2 December 2025
'''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytwalk

base_path = os.path.dirname(os.path.abspath(__file__))

# ============================== settings ==============================
var_scaling = 0.2   # variance scaling factor
a = 1.0             # alpha for dirichlet distribution (prior weights)
period = "hist"     # "hist", "MH" or "both"

# if "both", specify hist posterior (as prior for MH)
hist_params_path = os.path.join(base_path, "test_hist.txt")

iterations = 65e+5
thinning = 650
output_path = os.path.join(base_path, f"test_{period}")
# ======================================================================

# 1) LOAD DATA
# a) (PROXY) OBSERVATIONS
if period == "hist":
    proxy_file = os.path.join(base_path, "proxies_historical_1950-1999.csv")
elif period == "MH" or period == "both":
    proxy_file = os.path.join(base_path, "proxies_MH.csv")
else:
    raise ValueError("Invalid period. Use 'hist', 'MH', or 'both'.")

lnPIP_df = pd.read_csv(proxy_file)

lnPdIP = lnPIP_df["lnPdIP"].dropna().to_numpy()
lnPbIP = lnPIP_df["lnPbIP"].dropna().to_numpy()
proxy_data = np.append(lnPdIP, lnPbIP)

dino = lnPIP_df.dropna(subset=["lnPdIP"])
dino_indices = dino.index.to_numpy()
n_sites_dino = len(dino_indices)


# b) (MODEL) SIMULATIONS
if period == "hist":
    model_file = os.path.join(base_path, "models_historical_1950-1999.xlsx")
else:
    model_file = os.path.join(base_path, "models_MH.xlsx")

dino_df = pd.read_excel(model_file, sheet_name="dino")
bras_df = pd.read_excel(model_file, sheet_name="bras")

dino_model_data = dino_df.iloc[:, 1:].values
dino_model_data = dino_model_data[:, dino_indices]
bras_model_data = bras_df.iloc[:, 1:].values

model_data = np.concatenate((dino_model_data, bras_model_data), axis=1)
model_data = model_data.T

n_models = len(dino_df["model"])


# c) BAYSIC VARIANCE
dino_posterior_path = os.path.join(base_path, "server_dino_spavar_7900_2024-09-06_08-56-41.txt")
bras_posterior_path = os.path.join(base_path, "server_bras_spavar_7900_2024-09-06_08-56-41.txt")

dino_posterior = np.genfromtxt(dino_posterior_path, delimiter="\t")
dino_phi = dino_posterior[:, 2]
dino_phi_mean = np.mean(dino_phi)*var_scaling

bras_posterior = np.genfromtxt(bras_posterior_path, delimiter="\t")
bras_phi = bras_posterior[:, 2]
bras_phi_mean = np.mean(bras_phi)*var_scaling


# d) IF "BOTH", HISTORICAL PRIOR
if period == "both":
    hist_params = np.loadtxt(hist_params_path, delimiter="\t")
    hist_params = np.hstack((hist_params, (1 - np.sum(hist_params, axis=1)).reshape(-1, 1)))

    hist_params_mean = np.mean(hist_params, axis=0)
    hist_params_var = np.var(hist_params, axis=0)

    hist_ai = hist_params_mean**2 * (1-hist_params_mean) / hist_params_var - hist_params_mean
    a0 = np.sum(hist_ai)
    ai = a0*hist_params_mean



# 2) BAYESIAN WEIGHTING
def init():
    w = np.random.dirichlet(np.ones(n_models))
    w = w[:-1]
    return w


# Set alpha for dirichlet distribution (prior weights)
if period == "both":
    alpha = ai
else:
    alpha = np.full(n_models, a)

def prior(params):
    w_prior = np.sum((alpha-1)*np.log(params))
    return w_prior


# Set parameters for t distribution
t_a = 3
t_b = 4

def logLikelihood(params, proxy_data, model_data):
    weighted_sum = np.dot(model_data, params)
    res = -(2*t_a+1)/2 * np.sum(np.log(t_b + (proxy_data[:n_sites_dino]-weighted_sum[:n_sites_dino])**2 / (2*dino_phi_mean)))
    res += -(2*t_a+1)/2 * np.sum(np.log(t_b + (proxy_data[n_sites_dino:]-weighted_sum[n_sites_dino:])**2 / (2*bras_phi_mean)))
    return res


def obj(params):
    params = np.append(params, 1-params.sum())
    sum = logLikelihood(params, proxy_data, model_data) + prior(params)
    return -sum


def supp(params):
    if np.any(params < 0):
        return False
    elif np.sum(params) >= 1:
        return False
    else:
        return True


# Run t-walk
x0 = init()
xp0 = init()

twalk = pytwalk.pytwalk(n=len(init()), U=obj, Supp=supp)
twalk.Run(T=int(iterations), x0=x0, xp0=xp0, thin=thinning)

# Visualise convergence
plt.plot(twalk.Output[int(.5*len(twalk.Output)):, -1])
plt.title(f"iterations = {int(iterations)}, thinning = {thinning}")
plt.savefig(output_path + ".png")

# Save parameters (weights, excluding the last!)
parameters = twalk.Output[int(.5*len(twalk.Output)):, :-1]
df = pd.DataFrame(parameters)
df.to_csv(output_path + ".txt", sep="\t", index=False, header=False)
