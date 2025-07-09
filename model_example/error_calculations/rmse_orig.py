import numpy as np
import pandas as pd



df = pd.read_csv("all_results_5ns.csv")[:-3]
df.columns = ["mutant", "og_pred", "og_ci", "pred_std_error", "empirical", "quant_accuracy"]
n_bootstraps = 10000

emps = list(df["empirical"].astype("float64"))
og_preds = list(df["og_pred"])
og_cis = list(df["og_ci"])


print("Non bootstrapped original RMSE (5 ns, all systems): ")
print(np.sqrt(np.mean(np.array([(pred - emp) for pred, emp in zip(og_preds, emps)])**2)))
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


print("RMSE metrics")
print("Original calculations")
display_output(og_RMSEs)
print("------------------------------------------------------------")
print()


og_MAEs = bootstrap_mae(n_bootstraps, og_preds, og_cis, emps)
 

print()
print()
print("MAE metrics")
print("Original calculations")
display_output(og_MAEs)
print("------------------------------------------------------------")
print()

