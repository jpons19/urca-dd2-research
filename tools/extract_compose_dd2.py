#!/usr/bin/env python3
"""Extract the PT(GRDF2-DD2) inputs used by numerics/ from the original CompOSE archives.

Usage:
    python tools/extract_compose_dd2.py --cold cold.zip --general general.zip --out numerics/data/

Inputs (not redistributed): cold beta-equilibrated table, CompOSE EOS 302
(sha256 7e9d8acb11491ae64748be3f76b7cb6ad8d4724da2926ef505ea88d9de27f247, 73 922 bytes), and
general-purpose 3D table, CompOSE EOS 301
(sha256 163efbc9a59752c6d7e3e4257f34504c3df41af963c1d323384c8f2be2f02c3f, 190 451 886 bytes).

Outputs:
  DD2_cold_composition.csv     cold-table fractions, Fermi momenta and DU gaps (validation only);
  DD2_micro_first_two_T.npz    micro fields 10051:V_n, 10052:S_n, 11051:V_p, 11052:S_p on the
                               (nB, Yq) grid for the two lowest temperatures, order (T, Yq, nB, field).
"""
import argparse
import zipfile
from pathlib import Path
import numpy as np

HBARC = 197.3269804


def grid(z, name):
    x = z.read(name).split()
    n = int(x[1]); vals = np.array([float(y) for y in x[2:]], dtype=float)
    if len(vals) != n:
        raise ValueError(f'{name}: expected {n}, got {len(vals)}')
    return vals


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--cold', required=True, type=Path)
    a.add_argument('--general', required=True, type=Path)
    a.add_argument('--out', required=True, type=Path)
    args = a.parse_args(); args.out.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(args.cold) as z:
        nb = grid(z, 'eos.nb'); frac = np.zeros((len(nb), 4))
        with z.open('eos.compo') as f:
            for line in f:
                xs = line.split()
                if not xs:
                    continue
                j = int(xs[1]) - 1; count = int(xs[4]); pairs = xs[5:5 + 2 * count]
                sp = {int(pairs[2 * k]): float(pairs[2 * k + 1]) for k in range(count)}
                frac[j, :] = [sp.get(k, 0.) for k in [10, 11, 0, 1]]
        core = nb >= 0.16     # npe-mu neutrality only holds in the homogeneous core
        if np.max(np.abs(frac[core, 1] - frac[core, 2] - frac[core, 3])) > 1e-5:
            raise ValueError('core charge neutrality unexpectedly violated')
        p = HBARC * np.cbrt(3 * np.pi**2 * nb[:, None] * frac)
        dat = np.column_stack([nb, frac, p, p[:, 0] - p[:, 1] - p[:, 2], p[:, 0] - p[:, 1] - p[:, 3]])
        header = 'nB_fm3,Yn,Yp,Ye,Ymu,pFn_MeV,pFp_MeV,pFe_MeV,pFmu_MeV,DU_gap_e_MeV,DU_gap_mu_MeV'
        np.savetxt(args.out / 'DD2_cold_composition.csv', dat, delimiter=',', header=header, comments='')

    with zipfile.ZipFile(args.general) as z:
        nb = grid(z, 'eos.nb'); tg = grid(z, 'eos.t'); yq = grid(z, 'eos.yq')
        micro = np.full((2, len(yq), len(nb), 4), np.nan)
        with z.open('eos.micro') as f:
            for line in f:
                xs = line.split(); it, ib, iy, _ = map(int, xs[:4])
                if it > 2:
                    break
                if list(map(int, xs[4::2])) != [10051, 10052, 11051, 11052]:
                    raise ValueError('unexpected eos.micro fields')
                micro[it - 1, iy - 1, ib - 1] = np.array(xs[5::2], dtype=float)
        if not np.all(np.isfinite(micro)):
            raise ValueError('missing general-purpose micro entries')
        np.savez_compressed(args.out / 'DD2_micro_first_two_T.npz', nB=nb, Tgrid=tg[:2], Yqgrid=yq, fields=micro)
    print('Extracted DD2 cold composition and micro fields into', args.out)


if __name__ == '__main__':
    main()
