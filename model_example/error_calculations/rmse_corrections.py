import numpy as np
import pandas as pd


df = pd.read_csv("kelley_hu4d5_orig_corrected_5ns.csv")
df.columns = ["System", "og_pred", "og_ci", "corr_pred", "corr_ci", "std_error", "empirical", "r2"]

df_25ns = pd.read_csv("all_results_25ns.csv")
df_25ns.columns = ["system", "mutant", "og_pred", "og_ci", "std_err", "emp_value", "pred_err"]
df_25ns = df_25ns[df_25ns["system"] == "hu4D5-5"]

n_bootstraps = 10000

emps = list(df["empirical"].astype("float64"))
og_preds = list(df["og_pred"])
og_cis = list(df["og_ci"])
corr_preds = list(df["corr_pred"])
corr_cis = list(df["corr_ci"])


print("Non bootstrapped original RMSE (5 ns, only hu4D5-5 systems): ")
print(np.sqrt(np.mean(np.array([(pred - emp) for pred, emp in zip(og_preds, emps)])**2)))
print("Non bootstrapped corrected RMSE (only hu4D5-5 systems): ")
print(np.sqrt(np.mean(np.array([(pred - emp) for pred, emp in zip(corr_preds, emps)])**2)))
print("Non bootstrapped long RMSE (25 ns, no correction, only hu4D5-5 systems): ")
print(np.sqrt(np.mean(np.array([(pred - emp) for pred, emp in zip(df_25ns["og_pred"], df_25ns["emp_value"].astype("float64"))])**2)))
print()

def bootstrap_rmse(n_bootstraps, preds, spreads, emps, std_ci="ci"):
    RMSEs = []

    for i in range(n_bootstraps):
        boot_preds = []
        for center, spread in zip(preds, spreads):
            if (std_ci == "ci"):
                std = spread/1.96   
            else:
                std = spread
            boot_preds.append(np.random.normal(loc=center, scale=std, size=1)[0])

        errs = np.array([(pred - emp) for pred, emp in zip(boot_preds, emps)])
        rmse = np.sqrt(np.mean(errs**2))
        RMSEs.append(rmse)

    return RMSEs


def bootstrap_mae(n_bootstraps, preds, spreads, emps, std_ci="ci"):
    MAEs = []

    for i in range(n_bootstraps):
        boot_preds = []
        for center, spread in zip(preds, spreads):
            if (std_ci == "ci"):
                std = spread/1.96
            else:
                std = spread
            boot_preds.append(np.random.normal(loc=center, scale=std, size=1)[0])

        errs = np.array([(pred - emp) for pred, emp in zip(boot_preds, emps)])
        MAEs.append(np.mean(abs(errs)))

    return MAEs

def display_output(RMSE_arr):
    print("Summary stats")
    print(pd.Series(RMSE_arr).describe())

    rmse_cis = [np.percentile(RMSE_arr, 2.5), np.percentile(RMSE_arr, 97.5)]
    print(f"Low CI: {round(rmse_cis[0], 5)}")
    print(f"High CI: {round(rmse_cis[1], 5)}")

    ci_range = max(abs(np.mean(RMSE_arr) - rmse_cis))
    print(f"Total: {round(np.mean(RMSE_arr), 5)} +/- {round(ci_range, 5)}")

og_RMSEs = bootstrap_rmse(n_bootstraps, og_preds, og_cis, emps)
corr_RMSEs = bootstrap_rmse(n_bootstraps, corr_preds, corr_cis, emps)
long_RMSEs = bootstrap_rmse(
    n_bootstraps,
    list(df_25ns["og_pred"]),
    list(df_25ns["og_ci"]),
    list(df_25ns["emp_value"].astype("float64"))
    )


print("RMSE metrics")
print("Original calculations")
display_output(og_RMSEs)
print("------------------------------------------------------------")
print()
print("Corrected calculations")
display_output(corr_RMSEs)
print("------------------------------------------------------------")
print()
print("25 ns calculations")
display_output(long_RMSEs)


og_MAEs = bootstrap_mae(n_bootstraps, og_preds, og_cis, emps)
corr_MAEs = bootstrap_mae(n_bootstraps, corr_preds, corr_cis, emps)
long_MAEs = bootstrap_mae(
    n_bootstraps,
    list(df_25ns["og_pred"]),
    list(df_25ns["og_ci"]),
    list(df_25ns["emp_value"].astype("float64"))
    )
 

print()
print()
print("MAE metrics")
print("Original calculations")
display_output(og_MAEs)
print("------------------------------------------------------------")
print()
print("Corrected calculations")
display_output(corr_MAEs)
print("------------------------------------------------------------")
print()
print("25 ns calculations")
display_output(long_MAEs)

