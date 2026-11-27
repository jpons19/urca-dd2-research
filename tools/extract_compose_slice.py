#!/usr/bin/env python3
"""Extract the lowest-temperature slice of a CompOSE 3D general-purpose table.

Usage:
    python tools/extract_compose_slice.py --dir external_eos/BCP_NLRMF_N \
        --out numerics/data/BCP_NLRMF_N_T0.npz

Reads eos.nb, eos.yq, eos.t, eos.micro, eos.thermo and eos.compo (uncompressed,
CompOSE ASCII format) and keeps only temperature index 1.  Stored arrays:
  nB (n_nb), Yq (n_yq), T0 (MeV), masses (m_n, m_p) from eos.thermo,
  micro[n_yq, n_nb, k] for the micro indices listed in ``micro_idx``,
  mu_n, mu_p, mu_e [n_yq, n_nb] in MeV from eos.thermo,
  Yp [n_yq, n_nb] from eos.compo (index 11).
The original tables are not redistributed; only this derived slice is kept.
"""
import argparse
from pathlib import Path
import numpy as np


def grid(path):
    x = Path(path).read_text().split()
    n = int(x[1])
    vals = np.array([float(v) for v in x[2:2 + n]])
    if len(vals) != n:
        raise ValueError(f'{path}: expected {n} values, got {len(vals)}')
    return vals


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--dir', required=True, type=Path)
    a.add_argument('--out', required=True, type=Path)
    args = a.parse_args()
    d = args.dir
    nb, yq, tg = grid(d / 'eos.nb'), grid(d / 'eos.yq'), grid(d / 'eos.t')
    shape = (len(yq), len(nb))

    micro_idx, micro = None, None
    with open(d / 'eos.micro') as f:
        for line in f:
            x = line.split()
            it, ib, iy, n = map(int, x[:4])
            if it != 1:
                continue
            idx = [int(v) for v in x[4:4 + 2 * n:2]]
            if micro_idx is None:
                micro_idx = idx
                micro = np.full(shape + (len(idx),), np.nan)
            if idx != micro_idx:
                raise ValueError(f'micro index set changed at {(it, ib, iy)}: {idx}')
            micro[iy - 1, ib - 1] = [float(v) for v in x[5:5 + 2 * n:2]]

    with open(d / 'eos.thermo') as f:
        mn, mp = map(float, f.readline().split()[:2])
        mu_b = np.full(shape, np.nan); mu_q = np.full(shape, np.nan); mu_l = np.full(shape, np.nan)
        for line in f:
            x = line.split()
            it, ib, iy = map(int, x[:3])
            if it != 1:
                continue
            q3, q4, q5 = float(x[5]), float(x[6]), float(x[7])
            mu_b[iy - 1, ib - 1] = mn * (1 + q3)
            mu_q[iy - 1, ib - 1] = mn * q4
            mu_l[iy - 1, ib - 1] = mn * q5

    yp = np.full(shape, np.nan)
    with open(d / 'eos.compo') as f:
        for line in f:
            x = line.split()
            it, ib, iy = map(int, x[:3])
            if it != 1:
                continue
            n = int(x[4])
            sp = {int(x[5 + 2 * k]): float(x[6 + 2 * k]) for k in range(n)}
            yp[iy - 1, ib - 1] = sp.get(11, 0.0)

    for name, arr in (('micro', micro), ('mu_b', mu_b), ('Yp', yp)):
        if not np.all(np.isfinite(arr)):
            raise ValueError(f'missing entries in {name} for T index 1')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, nB=nb, Yq=yq, T0=tg[0], masses=np.array([mn, mp]),
                        micro_idx=np.array(micro_idx), micro=micro,
                        mu_n=mu_b, mu_p=mu_b + mu_q, mu_e=mu_l - mu_q, Yp=yp)
    print(f'wrote {args.out}: nB {nb[0]:.3g}-{nb[-1]:.3g} ({len(nb)}), Yq {yq[0]}-{yq[-1]} ({len(yq)}), '
          f'T0={tg[0]} MeV, micro {micro_idx}')


if __name__ == '__main__':
    main()
