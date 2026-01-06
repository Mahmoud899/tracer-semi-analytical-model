#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from simulation_options import simulateRoseNDSSinglePoro
from sampling_algorithms import lhs_sample
from metrics_options import least_squares_error, calcR2
from Rose_data import nds_true_conc

# --- choose how many CPUs to use ---
MAX_WORKERS = 120  # set this to the number of processes

# --- worker: runs one simulation and returns results ---
def _run_one(it, mean_residence_time_1, peclet_number_1, fractional_recovery_1,
                 mean_residence_time_2, peclet_number_2):
    prediction_1 = simulateRoseNDSSinglePoro(
        mean_residence_time=mean_residence_time_1,
        peclet_number=peclet_number_1
    )
    prediction_2 = simulateRoseNDSSinglePoro(
        mean_residence_time=mean_residence_time_2,
        peclet_number=peclet_number_2
    )
    
    fractional_recovery_2 = 1 - fractional_recovery_1

    prediction = prediction_1 * fractional_recovery_1 + prediction_2 * fractional_recovery_2
    mse = least_squares_error(prediction, nds_true_conc, mean=True)
    r2_value = calcR2(prediction, nds_true_conc)
    return it, mse, r2_value

def main():
    # Sample the simulation experiments
    parameters = {
        'mean_residence_time_1': [1, 50],
        'peclet_number_1': [1, 100],
        'mean_residence_time_2': [1, 50],
        'peclet_number_2': [1, 100],
        'fractional_recovery': [0.01, 1]
    }
    sample = lhs_sample(parameters, 1024*8)

    n = len(sample['iteration'])
    data_placeholder = sample.copy()

    # Pre-allocate results by iteration index
    mse_arr = np.empty(n, dtype=float)
    r2_arr  = np.empty(n, dtype=float)

    # Build task list (one per iteration)
    tasks = []
    for it in data_placeholder['iteration']:
        mrt_1  = data_placeholder['mean_residence_time_1'][it]
        pec_1  = data_placeholder['peclet_number_1'][it]
        mrt_2  = data_placeholder['mean_residence_time_2'][it]
        pec_2  = data_placeholder['peclet_number_2'][it]
        frec = data_placeholder['fractional_recovery'][it]
        tasks.append((it, mrt_1, pec_1, mrt_2, pec_2, frec))

    # Run in parallel with a single progress bar
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as exe:
        futures = [exe.submit(_run_one, it, mrt_1, pec_1, frec, mrt_2, pec_2) for (it, mrt_1, pec_1, mrt_2, pec_2, frec) in tasks]
        with tqdm(total=len(futures), desc='Simulating iterations ..', unit='iter', position=0, leave=True) as pbar:
            for fut in as_completed(futures):
                it, mse, r2_value = fut.result()
                mse_arr[it] = mse
                r2_arr[it]  = r2_value
                pbar.update(1)

    # Assemble DataFrame and save
    data_placeholder['MSE'] = mse_arr.tolist()
    data_placeholder['R2']  = r2_arr.tolist()

    df = pd.DataFrame(data_placeholder)
    df.to_csv('SinglePorosityTwoFractureSimulationRuns.csv', index=False)

if __name__ == '__main__':
    main()
