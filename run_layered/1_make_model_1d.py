#!/usr/bin/env python3
import numpy as np
from typing import List, Optional

def find_index(x, val):
    """
    Find index k such that x[k] <= val < x[k+1].
    If val < x[0] returns 0.
    If val >= x[-1] returns len(x)-2 (so k+1 is last index).
    This mirrors the behavior in the C binary-search code.
    x is assumed ascending and len(x) >= 2.
    """
    if val < x[0]:
        return 0
    # searchsorted returns first index i where x[i] > val (if side='right')
    i = np.searchsorted(x, val, side='right')
    k = max(0, i - 1)
    # ensure k not beyond len(x)-2
    if k >= len(x) - 1:
        return len(x) - 2
    return k
import math

def create_nugrid(n: int, length: float, dx: float, x: Optional[List[float]] = None) -> float:
    """
    Generate nonuniform grid using geometric progression.
    Determine the optimal ratio r by root finding using fixed point iteration.
    
    Parameters:
    n: int - number of intervals (n+1 points/nodes)
    length: float - total length of the grid
    dx: float - first interval size
    
    Returns:
    float - the ratio r for geometric progression
    """
    
    eps = 1e-15
    
    # If uniform grid is sufficient
    if abs(n * dx - length) < eps:
        return 1.0
    
    # Fixed point iteration to find the root "r" of the equation:
    # length = dx * (r^n - 1) / (r - 1)
    r = 1.1
    rr = 1.0
    
    while True:
        rr = (length * (r - 1.0) / dx + 1.0) ** (1.0 / n)
        if abs(rr - r) < eps:
            break
        r = rr
    x[0] = 0
    for i in range(1, n + 1):
        x[i] = (math.pow(r, i) - 1.0) * dx / (r - 1.0)
    
    return r

def main(argv=None):
    n = 50
    x = np.linspace(0, 10000, n+1)
    #r = create_nugrid(n, 10000, 100, x)
    x1_vec = np.r_[-x[::-1], x[1:]]
    #x1_vec = np.arange(-10e3,10e3,200)
    x2_vec = x1_vec
    x3_vec = np.linspace(0, 5000, 101)
    print(x3_vec)

    depth = np.array([1000, 2300, 2400])
    resis = np.array([0.3, 1, 50, 2])
    aniso = np.array([1, 1, 1.414, 1])
    nanomaly = 0

    # midpoints (cell centers)
    x1s = 0.5 * (x1_vec[:-1] + x1_vec[1:])
    x2s = 0.5 * (x2_vec[:-1] + x2_vec[1:])
    x3s = 0.5 * (x3_vec[:-1] + x3_vec[1:])

    ndepth = depth.size
    n1 = x1_vec.size-1
    n2 = x2_vec.size-1
    n3 = x3_vec.size-1
    
    # allocate rho arrays with shape (n3, n2, n1)
    rho11 = np.full((n3, n2, n1), resis[0], dtype=np.float32)
    rho22 = np.full((n3, n2, n1), resis[0], dtype=np.float32)
    # rho33 scaled by aniso^2 (top layer uses aniso[0])
    rho33 = np.full((n3, n2, n1), resis[0] * (aniso[0] ** 2), dtype=np.float32)

    # apply depth layers
    # For each interface depth[j], find i3min such that x3s[i3min] <= depth[j] < x3s[i3min+1]
    for j in range(ndepth):
        i3min = find_index(x3s, depth[j])
        # in C they do for(i3=i3min+1; i3<mod->n3; i3++) => layers below interface are set
        for i3 in range(i3min + 1, n3):
            rho11[i3, :, :] = resis[j+1]
            rho22[i3, :, :] = resis[j+1]
            rho33[i3, :, :] = resis[j+1] * (aniso[j+1] ** 2)

    # apply anomalies
    if nanomaly > 0:
        for j in range(nanomaly):
            x1min_j, x1max_j = x1min[j], x1max[j]
            x2min_j, x2max_j = x2min[j], x2max[j]
            x3min_j, x3max_j = x3min[j], x3max[j]
            # create boolean masks for each axis
            mask_x1 = (x1s >= x1min_j) & (x1s <= x1max_j)
            mask_x2 = (x2s >= x2min_j) & (x2s <= x2max_j)
            mask_x3 = (x3s >= x3min_j) & (x3s <= x3max_j)
            # iterate over logical indices - broadcasting to 3D
            # we want to set rho[..., i2, i1] where all three masks true
            # Build 3D mask with broadcasting
            mask3d = mask_x3[:, None, None] & mask_x2[None, :, None] & mask_x1[None, None, :]
            rho11[mask3d] = anomaly[j]
            rho22[mask3d] = anomaly[j]
            rho33[mask3d] = anomaly[j]

    # write binary files (float32), match C's layout:
    # fx1 (n1+1 floats), fx2 (n2+1), fx3 (n3+1)
    x1_vec.astype(np.float32).tofile("fx1")
    x2_vec.astype(np.float32).tofile("fx2")
    x3_vec.astype(np.float32).tofile("fx3")

    # frho... write as contiguous float32 array. The C code likely writes with i1 fastest
    # ndarray with shape (n3, n2, n1) and C-order has i1 the fastest index - so tofile is compatible.
    rho11.astype(np.float32).tofile("frho11")
    rho22.astype(np.float32).tofile("frho22")
    rho33.astype(np.float32).tofile("frho33")

    print("Wrote fx1, fx2, fx3, frho11, frho22, frho33")

if __name__ == "__main__":
    main()
