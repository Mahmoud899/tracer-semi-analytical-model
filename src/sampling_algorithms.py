import os
import numpy as np
from typing import Dict, List, Tuple

def lhs_sample(params: Dict[str, Tuple[float, float]], n: int, seed: int = None) -> Dict[str, List[float]]:
    """
    Latin Hypercube Sampling (uniform within each stratum) over the ranges in `params`.

    params: {"param_name": [min, max], ...}
    n     : number of samples (strata per dimension)
    seed  : optional RNG seed for reproducibility

    returns: {"iteration": [0..n-1], "param1": [...], "param2": [...], ...}
    """
    if n <= 0:
        raise ValueError("n must be a positive integer.")
    rng = np.random.default_rng(seed)

    names = list(params.keys())
    d = len(names)
    samples = np.empty((n, d), dtype=float)

    # Build LHS in [0,1]^d
    for j, name in enumerate(names):
        lo, hi = params[name]
        if hi < lo:
            raise ValueError(f"Range for '{name}' must satisfy min <= max.")
        # permutation of strata 0..n-1
        strata = rng.permutation(n)
        # one point uniformly within each stratum
        u = rng.random(n)
        x = (strata + u) / n  # in (0,1)
        # scale to [lo, hi]
        samples[:, j] = lo + x * (hi - lo)

    # Package as dict of lists
    out = {"iteration": list(range(n))}
    for j, name in enumerate(names):
        out[name] = samples[:, j].tolist()
    return out



