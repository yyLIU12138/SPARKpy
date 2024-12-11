# Created based on SPARK package (Version: 1.0.2) implemented in R
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.genmod.families import Poisson
from statsmodels.tools import add_constant
from scipy.sparse import csr_matrix, csc_matrix
from scipy.stats import chi2
from scipy.stats import cauchy
from scipy.linalg import eigh
from statsmodels.stats.multitest import multipletests

from .utils import cal_kernel_params, gaussian_kernel, cosine_kernel


class SPARKpy:

    def __init__(self, adata, use_rep='spatial', cov_mtx=None, max_iter=500, tol=1e-5, fit='Poisson', verbose=True):

        self.adata = adata  # data to be analyzed
        self.use_rep = use_rep  # key of representation to use (e.g. spatial coordinates)
        self.cov_mtx = cov_mtx  # adjusted covariate matrix (batch indices etc.)
        if self.cov_mtx is not None:
            if self.cov_mtx.ndim == 1:
                self.cov_mtx = self.cov_mtx.reshape(-1, 1)
            else:
                self.cov_mtx = self.cov_mtx
        self.max_iter = max_iter  # max interation for updating
        self.tol = tol  # threshold for convergence
        self.fit = fit  # model type

        self.verbose = verbose

        # calculate library size
        if isinstance(adata.X, (csr_matrix, csc_matrix)):
            self.lib_size = np.sum(adata.X.toarray(), axis=1)
        else:
            self.lib_size = np.sum(adata.X, axis=1)
        self.loglib = np.log(self.lib_size)

        self.n_spots = self.adata.shape[0]  # number of spots
        self.n_genes = self.adata.shape[1]  # number of genes
        if self.cov_mtx is not None:
            self.n_cov = self.cov_mtx.shape[1]  # number of adjusted covariates
        else:
            self.n_cov = 0

        if self.verbose:
            print("## ===== SPARKpy INPUT INFORMATION ==== \n")
            print(f"## number of total samples: {self.n_spots}\n")
            print(f"## number of total features: {self.n_genes}\n")
            print(f"## number of adjusted covariates: {self.n_cov}\n")

    def fit_spark(self):

        if self.fit.lower() == 'poisson':

            if self.verbose:
                print("## ===== fitting Poisson model ==== \n")

            valid_genes = []
            p_all_list = []
            p_T_list = []

            mtx = np.array(self.adata.obsm[self.use_rep])
            diff = mtx[:, np.newaxis, :] - mtx[np.newaxis, :, :]
            dist_mtx = np.sqrt(np.sum(diff ** 2, axis=2))
            param_list = cal_kernel_params(dist_mtx=dist_mtx, min_idx=3, counts=5)

            if self.verbose:
                print(f"## ===== preparing kernels ==== \n")
            kernel_list = []
            for param in param_list:
                kernel_g = gaussian_kernel(dist_mtx=dist_mtx, sigma=param)
                kernel_list.append(kernel_g)
                kernel_c = cosine_kernel(dist_mtx=dist_mtx, phi=param)
                kernel_list.append(kernel_c)

            for k, kernel in enumerate(kernel_list):
                eigval, eigvec = eigh(kernel)
                eigval[eigval < 1e-8] = 1e-8
                kernel_list[k] = np.dot(eigvec * eigval, eigvec.T)

            for g in range(self.n_genes):

                if self.verbose:
                    print(f"## ===== fitting model for the {g}-th gene ==== \n")
                try:
                    y_tilde, P = self.fit_poisson(g)

                    if np.isnan(y_tilde).any() or np.isnan(P).any():
                        # if self.verbose:
                        print(f"Skipping gene {g} due to NaN values.")
                        continue
                    valid_genes.append(list(self.adata.var_names)[g])

                except Exception as e:
                    # if self.verbose:
                    print(f"Error in fitting model for gene {g}: {e}. Skipping this gene.")
                    continue

                if self.verbose:
                    print(f"## ===== testing the {g}-th gene ==== \n")

                if self.verbose:
                    print(f"## ===== calculating p-values ==== \n")
                p_value_list = []
                for kernel in kernel_list:
                    p_value = self.test_poisson(kernel, y_tilde, P)
                    p_value_list.append(p_value)

                if self.verbose:
                    print(f"## ===== combining p-values using Cauchy combination rule ==== \n")
                p_T = self.cauchy_combination(p_value_list)

                p_all_list.append(p_value_list)
                p_T_list.append(p_T)

        else:

            raise ValueError("Only Poisson version of SPARK is implemented in current version !")

        if self.verbose:
            print(f"## ===== calculating adjusted p-values according to B-Y procedure ==== \n")

        # adj_p_values = multipletests(np.array(p_T_list), method='fdr_by')[1]

        cm = np.sum(1 / np.arange(1, len(p_T_list) + 1))
        sorted_indices = np.argsort(p_T_list)
        adj_p_values = np.array(p_T_list) * len(p_T_list) * cm / (sorted_indices + 1)
        adj_p_values = list(np.clip(adj_p_values, 0, 1.0))

        head_list = []
        p_all_array = np.array(p_all_list)
        for i in range(len(param_list)):
            head_list.append(f'gau{i}')
            head_list.append(f'cos{i}')

        p_value_df = pd.DataFrame(p_all_array, columns=head_list, index=valid_genes)
        p_value_df['p_values'] = p_T_list
        p_value_df['adj_p_values'] = adj_p_values

        return p_value_df

    def fit_poisson(self, g):

        # covariate matrix
        if not self.n_cov:
            X = np.ones((self.n_spots, 1))  # intercept only
        else:
            X = add_constant(self.cov_mtx)

        y = self.adata.X[:, g].toarray().squeeze()  # counts

        # fit GLM with Poisson family (link function: log, offset: log library size)
        model0 = sm.GLM(y, X, family=Poisson(link=sm.families.links.Log()),
                        offset=self.loglib).fit()

        tau2 = 1.  # tau1: spatial related / tau2: random
        beta = model0.params  # * k
        mu = model0.fittedvalues  # * n
        eps = np.zeros(self.n_spots)  # * n

        H_diag = 1 / mu + tau2  # * n
        H_inv_diag = 1 / H_diag  # * n
        H_inv = np.diag(H_inv_diag)  # n * n

        t = 0
        while(t < self.max_iter):
            t += 1
            step_size = 1

            if X.shape[1] > 1:

                # pseudo data
                y_tilde = (X @ beta.reshape(-1, 1)).squeeze() + eps + (y - mu) / mu  # * n
                # y_tilde[np.abs(y_tilde) > 1e3] = np.median(y_tilde)

                P = H_inv - H_inv @ X @ np.linalg.inv(X.T @ H_inv @ X) @ X.T @ H_inv  # n * n

                Py = P @ y_tilde.reshape(-1, 1)  # n * 1
                PPy = P @ Py  # n * 1
                PPPy = P @ PPy  # n * 1

                tau2_new = max(0, (tau2 + step_size * 1 / (y_tilde.reshape(1, -1) @ PPPy).squeeze() *
                                   ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P))))
                z = 0
                while tau2_new == 0 and z < self.max_iter:
                    z += 1
                    step_size *= 0.5
                    tau2_new = max(0, (tau2 + step_size * 1 / (y_tilde.reshape(1, -1) @ PPPy).squeeze() *
                                       ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P))))
                # tau2_new = max(0, (tau2 + 1 / (y_tilde.reshape(1, -1) @ PPPy).squeeze() *
                #             ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P)) / self.n_spots))  # scalar
                # tau2_new = max(0, (tau2 + tau2 ** 2 *
                #             ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P)) / self.n_spots)) # same as the modified version of SPARK 2019

                H_diag = 1 / mu + tau2_new  # * n
                H_inv_diag = 1 / H_diag  # * n
                H_inv = np.diag(H_inv_diag)  # n * n

                # beta_new = np.linalg.inv(X.T @ H_inv @ X) @ X.T @ H_inv @ y_tilde  # * k
                beta_new = (np.linalg.inv(np.sum(np.array([H_inv[i, i] * np.outer(X[i], X[i]) for i in range(X.shape[0])]), axis=0)) @
                            np.sum((X * (H_inv_diag * y_tilde)[:, np.newaxis]).T, axis=1))

                eps_new = tau2_new * H_inv_diag * (y_tilde - (X @ beta_new.reshape(-1, 1)).squeeze())  # * n

                mu = self.lib_size * np.exp((X @ beta_new.reshape(-1, 1)).squeeze() + eps_new)  # * n

                if (np.max(np.abs(beta_new - beta) / (np.abs(beta) + np.abs(beta_new) + self.tol)) < self.tol and
                        np.max(np.abs(eps_new - eps) / (np.abs(eps) + np.abs(eps_new) + self.tol)) < self.tol):
                # if (np.max(np.abs(beta_new - beta)) < self.tol and
                #         np.max(np.abs(eps_new - eps) / (np.abs(eps) + np.abs(eps_new) + self.tol)) < self.tol):
                    tau2 = tau2_new
                    beta = beta_new
                    eps = eps_new
                    y_tilde = (X @ beta.reshape(-1, 1)).squeeze() + eps + (y - mu) / mu  # * n
                    # y_tilde[np.abs(y_tilde) > 1e3] = np.median(y_tilde)
                    P = H_inv - H_inv @ X @ np.linalg.inv(X.T @ H_inv @ X) @ X.T @ H_inv  # n * n
                    break

                tau2 = tau2_new
                beta = beta_new
                eps = eps_new
                if t == self.max_iter:
                    y_tilde = (X @ beta.reshape(-1, 1)).squeeze() + eps + (y - mu) / mu  # * n
                    # y_tilde[np.abs(y_tilde) > 1e3] = np.median(y_tilde)
                    P = H_inv - H_inv @ X @ np.linalg.inv(X.T @ H_inv @ X) @ X.T @ H_inv  # n * n
                    # print(f'{g}-th gene unconverged!')

            elif X.shape[1] == 1:

                # pseudo data
                y_tilde = (X * beta).squeeze() + eps + (y - mu) / mu  # * n
                # y_tilde[np.abs(y_tilde) > 1e3] = np.median(y_tilde)

                P = np.diag(H_inv_diag) - 1 / np.sum(H_inv_diag) * np.outer(H_inv_diag, H_inv_diag)  # n * n

                Py = P @ y_tilde.reshape(-1, 1)  # n * 1
                PPy = P @ Py  # n * 1
                PPPy = P @ PPy  # n * 1

                tau2_new = max(0, (tau2 + step_size * 1 / (y_tilde.reshape(1, -1) @ PPPy).squeeze() *
                                   ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P))))
                z = 0
                while tau2_new == 0 and z < self.max_iter:
                    z += 1
                    step_size *= 0.5
                    tau2_new = max(0, (tau2 + step_size * 1 / (y_tilde.reshape(1, -1) @ PPPy).squeeze() *
                                       ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P))))
                # tau2_new = max(0, (tau2 + 1 / (y_tilde.reshape(1, -1) @ PPPy).squeeze() *
                #             ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P)) / self.n_spots))  # scalar
                # tau2_new = max(0, (tau2 + tau2 ** 2 *
                #             ((y_tilde.reshape(1, -1) @ PPy).squeeze() - np.trace(P)) / self.n_spots))

                H_diag = 1 / mu + tau2_new  # * n
                H_inv_diag = 1 / H_diag  # * n
                # H_inv = np.diag(H_inv_diag)  # n * n

                beta_new = 1 / np.sum(H_inv_diag) * np.sum(H_inv_diag * y_tilde)  # scalar

                eps_new = tau2_new * H_inv_diag * (y_tilde - (X * beta_new).squeeze())  # * n

                mu = self.lib_size * np.exp((X * beta_new).squeeze() + eps_new)  # * n

                if (np.max(np.abs(beta_new - beta) / (np.abs(beta) + np.abs(beta_new) + self.tol)) < self.tol and
                        np.max(np.abs(eps_new - eps) / (np.abs(eps) + np.abs(eps_new) + self.tol)) < self.tol):
                # if (np.max(np.abs(beta_new - beta)) < self.tol and
                #         np.max(np.abs(eps_new - eps) / (np.abs(eps) + np.abs(eps_new) + self.tol)) < self.tol):
                    tau2 = tau2_new
                    beta = beta_new
                    eps = eps_new
                    y_tilde = (X * beta).squeeze() + eps + (y - mu) / mu  # * n
                    # y_tilde[np.abs(y_tilde) > 1e3] = np.median(y_tilde)
                    P = np.diag(H_inv_diag) - 1 / np.sum(H_inv_diag) * np.outer(H_inv_diag, H_inv_diag)  # n * n
                    break

                tau2 = tau2_new
                beta = beta_new
                eps = eps_new
                if t == self.max_iter:
                    y_tilde = (X * beta).squeeze() + eps + (y - mu) / mu  # * n
                    # y_tilde[np.abs(y_tilde) > 1e3] = np.median(y_tilde)
                    P = np.diag(H_inv_diag) - 1 / np.sum(H_inv_diag) * np.outer(H_inv_diag, H_inv_diag)  # n * n
                    # print(f'{g}-th gene unconverged!')

        return y_tilde, P

    def test_poisson(self, K, y_tilde, P):

        PK = P @ K

        trace_PKP = np.sum(PK * P)  # not included in the method, but included in the R version of SPARK

        es = np.trace(PK) / 2  # mean for S0
        I = np.trace(PK @ PK) / 2  # variance for S0
        I -= 0.5 * trace_PKP ** 2 / np.sum(P * P)  # not included in the method, but included in the R version of SPARK

        k_hat = 0.5 * I / es
        v_hat = 2 * es ** 2 / I

        S0 = (1/2 * y_tilde.reshape(1, -1) @ PK @ P @ y_tilde.reshape(-1, 1)).squeeze()
        S = S0 / k_hat

        p_value = 1 - chi2.cdf(S, df=v_hat)
        if np.isnan(p_value):
            print("Error encountered when calculating p-values !")

        return p_value

    def cauchy_combination(self, p_value_list, weight_list=None):
        # print(p_value_list)

        # avoid extreme values, thresholds are set the same as SPARK
        p_value_list = [(5.55e-17 if p == 0 else p) for p in p_value_list]
        p_value_list = [(0.99 if (1 - p) < 1e-3 else p) for p in p_value_list]

        if weight_list:
            weights = np.array(weight_list)
        else:
            weights = np.ones(len(p_value_list)) / len(p_value_list)
        T = np.sum(weights * np.tan((0.5 - np.array(p_value_list)) * np.pi))
        p_T = 1 / 2 - np.arctan(T / np.sum(weights)) / np.pi
        # p_T = 1.0 - cauchy.cdf(T)
        # print(p_T)
        if p_T <= 0:
            p_T = 5.55e-17
        if np.isnan(p_T):
            print("Error encountered when calculating p-values !")

        return p_T


















