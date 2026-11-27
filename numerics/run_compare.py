#!/usr/bin/env python3
"""Multi-EOS comparison of Urca emissivities versus baryon density.

From repository root:
    python -m numerics.run_compare [--points 120] [--eos DD2 BCP BSk22]

Writes results/<EOS>_urca_T9.csv, the comparison figure and an EOS-inputs figure.
Uses the derived inputs in numerics/data/; the TOV reference densities are cached
in numerics/data/stars.json (recomputing them needs the CompOSE tables in external_eos/).
"""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from . import eos_inputs as E
from . import rates as R
from . import tov
from . import mmu as M
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
RANGES = {'DD2': (0.12, 0.98), 'BCP': (0.12, 1.10), 'BSk22': (0.12, 1.20)}
LABELS = {'DD2': 'PT(GRDF2-DD2) [RMF]', 'BCP': 'BCP(NLRMF-N) [RMF]', 'BSk22': 'PCP(BSk22) [Skyrme]'}
STYLE = {'FM': ('Friman–Maxwell / NSCool', '--', 'C0'),
         'Bottaro': ('Bottaro triangular', '-', 'C1'),
         'Bottaro_nontri': ('Bottaro nontriangular', '-.', 'C2'),
         'Shternin': ('Shternin 0.6R_A+0.2R_B', '-', 'C3'),
         'NWA': ('Alford NWA (adapted)', '-', 'C4'),
         'Sedrakian': ('Sedrakian proton width', '-', 'C5'),
         'MMU': ('Voskresensky MMU (ratio form, $n_c^\\pi$=3$n_0$)', '-', 'C6'),
         'MMU_direct': ('MMU direct form (SVSWW97 Eq. 1)', ':', 'C6'),
         'MMU_unsoftened': ('MMU pion-mediated term, vacuum-like gap (ours)', '-.', 'C7'),
         'DU': ('Direct Urca (sharp)', ':', 'k')}
TOV_INPUT = {'DD2': ('DD2_cold', 'eos.thermo', 'eos.nb'), 'BCP': ('BCP_NLRMF_N', 'eos.thermo.ns', 'eos.nb.ns'),
             'BSk22': ('BSk22', 'eos.thermo', 'eos.nb')}
STATE_KEYS = ['nb', 'Yn', 'Yp', 'Ye', 'Ym', 'pn', 'pp', 'pe', 'pm', 'mue', 'mLn', 'mLp', 'mnD', 'mpD',
              'gap_e', 'gap_mu', 'beta_residual', 'charge_residual']


def thresholds(rows, key):
    """Densities where a DU momentum gap changes sign (linear interpolation)."""
    x = np.array([r['nb'] for r in rows]); g = np.array([r[key] for r in rows], dtype=float)
    out = []
    for i in range(len(x) - 1):
        if np.isfinite(g[i]) and np.isfinite(g[i + 1]) and g[i] > 0 >= g[i + 1]:
            out.append(x[i] + (x[i + 1] - x[i]) * g[i] / (g[i] - g[i + 1]))
    return out


def stars(name):
    """Central densities of 1.4 Msun and Mmax stars, and DU-onset masses (cached)."""
    cache = E.DATA / 'stars.json'
    data = json.loads(cache.read_text()) if cache.exists() else {}
    if name not in data:
        sub, th, nbf = TOV_INPUT[name]
        r = tov.reference_densities(E.EXTERNAL / sub, th, nbf)
        seq = r.pop('sequence')
        r['sequence_nc_R_M'] = seq.tolist()
        data[name] = r
        OUT.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(data, indent=1))
    return data[name]


def mass_at(star, nc):
    s = np.array(star['sequence_nc_R_M']); i = int(np.argmax(s[:, 2]))
    return float(np.interp(nc, s[:i + 1, 0], s[:i + 1, 2])) if nc <= s[i, 0] else np.nan


def scan(name, npoints):
    table = E.LOADERS[name]()
    lo, hi = RANGES[name]
    lo, hi = max(lo, table.nb_range[0]), min(hi, table.nb_range[1])
    rows = []
    for nb in np.linspace(lo, hi, npoints):
        d = table.state(float(nb))
        r = {k: d[k] for k in STATE_KEYS}
        for model, (qe, qm) in R.all_rates(d).items():
            r[f'{model}_e'], r[f'{model}_mu'] = qe, qm
            r[f'{model}_total'] = qe + qm
        for model, q in M.all_mmu(d).items():          # totals only (e + mu)
            r[f'{model}_total'] = q
        rows.append(r)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f'{name}_urca_T9.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(rows)
    return rows


def plot(results):
    names = list(results)
    fig, axes = plt.subplots(2, len(names), figsize=(5.2 * len(names), 10.2), sharex='col',
                             gridspec_kw={'height_ratios': [2.2, 1], 'hspace': .06, 'wspace': .22})
    axes = np.atleast_2d(axes).reshape(2, len(names))
    for j, name in enumerate(names):
        rows = results[name]; x = np.array([r['nb'] for r in rows])
        ref = np.array([r['FM_total'] for r in rows])
        ax, ar = axes[0, j], axes[1, j]
        for model, (label, ls, c) in STYLE.items():
            y = np.array([r[f'{model}_total'] for r in rows], dtype=float)
            if not np.any(np.isfinite(y) & (y > 0)):
                continue
            yp = np.where(y > 0, y, np.nan)
            ax.semilogy(x, yp, ls=ls, color=c, lw=2, label=label)
            if model != 'FM':
                ar.semilogy(x, yp / ref, ls=ls, color=c, lw=1.8)
            if model == 'MMU':
                y2 = np.array([r['MMU_nc2_total'] for r in rows], dtype=float)
                ax.fill_between(x, yp, y2, color=c, alpha=.15, lw=0, label=r'MMU: $n_c^\pi$ 3$n_0$ → 2$n_0$')
                ar.fill_between(x, yp / ref, y2 / ref, color=c, alpha=.15, lw=0)
            if model in ('Bottaro', 'Bottaro_nontri'):
                y1 = np.array([r[f'{model}_noSRC_total'] for r in rows], dtype=float)
                ax.fill_between(x, yp, y1, color=c, alpha=.18, lw=0,
                                label=r'$\beta_n$: 0.68 → 1' if model == 'Bottaro' else None)
                ar.fill_between(x, yp / ref, y1 / ref, color=c, alpha=.18, lw=0)
        for xe in thresholds(rows, 'gap_e'):
            for a in (ax, ar):
                a.axvline(xe, color='0.35', lw=1, ls=(0, (4, 3)))
        for xm in thresholds(rows, 'gap_mu'):
            for a in (ax, ar):
                a.axvline(xm, color='0.6', lw=1, ls=(0, (1, 2)))
        st = stars(name)
        for key, lab in (('nc_1.4', r'$n_c(1.4M_\odot)$'), ('nc_Mmax', r'$n_c(M_{\max})$')):
            if key in st:
                ax.axvline(st[key], color='tab:brown', lw=1.1, alpha=.7)
                ax.text(st[key], 3e27, ' ' + lab, rotation=90, va='top', fontsize=8, color='tab:brown')
        if st.get('nc_Mmax', 9) < x[-1]:
            for a in (ax, ar):
                a.axvspan(st['nc_Mmax'], x[-1], color='0.85', alpha=.5, lw=0)
        ax.set_title(LABELS[name] + '\n' + rf"$M_{{\max}}$={st['Mmax']:.2f} $M_\odot$", fontsize=11)
        ax.set_ylim(1e17, 1e28); ar.set_ylim(1e-4, 1e8)
        ar.axhline(1, color='0.3', lw=.8, ls=':')
        ar.set_xlabel(r'$n_B$ (fm$^{-3}$)')
        for a in (ax, ar):
            a.grid(alpha=.22, which='both')
        if j == 0:
            ax.set_ylabel(r'$Q_\nu$ at $T=10^9$ K (erg cm$^{-3}$ s$^{-1}$)')
            ar.set_ylabel('ratio to FM')
            pass
    fig.text(.06, .008, 'Neutron-branch MU only (FM, Bottaro, Shternin); e + μ channels summed. Dashed grey: electron DU threshold; '
             'dotted grey: muon DU threshold.\nShternin undefined at/above threshold; adapted NWA needs a Dirac mass (RMF only); Sedrakian: nonrelativistic proton spectral width, μ channel added by us. MMU: published pion-gap/vertex fits (n0=0.16), not EOS-consistent.\n'
             'Shaded: Bottaro without the inherited β=0.68 (possible overlap with ρ exchange). Differences between prescriptions are not an error band. Brown lines: central densities from TOV; grey: beyond $M_{\\max}$.', fontsize=8.5, color='.3')
    seen = {}
    for a in axes[0]:
        for h, l in zip(*a.get_legend_handles_labels()):
            seen.setdefault(l, h)
    handles, labels = list(seen.values()), list(seen)
    fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, 0.085))
    fig.subplots_adjust(left=.07, right=.985, top=.94, bottom=.25)
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'urca_multi_eos_T9.{ext}', dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10), sharex=True)
    for j, name in enumerate(names):
        rows = results[name]; x = np.array([r['nb'] for r in rows]); c = f'C{j}'
        axes[0].plot(x, [r['Yp'] for r in rows], color=c, lw=2, label=f'{LABELS[name]} $Y_p$')
        axes[0].plot(x, [r['Ym'] for r in rows], color=c, lw=1.2, ls='--', label=r'$Y_\mu$')
        axes[1].plot(x, [r['mLn'] for r in rows], color=c, lw=2, label=f'{name} $m_{{L,n}}$')
        axes[1].plot(x, [r['mLp'] for r in rows], color=c, lw=1.2, ls='--', label=f'{name} $m_{{L,p}}$')
        if np.any(np.isfinite([r['mnD'] for r in rows])):
            axes[1].plot(x, [r['mnD'] for r in rows], color=c, lw=1, ls=':', label=f'{name} $m_D$')
        axes[2].plot(x, [r['gap_e'] for r in rows], color=c, lw=2, label=f'{name} e')
        axes[2].plot(x, [r['gap_mu'] for r in rows], color=c, lw=1.2, ls='--', label=f'{name} μ')
    axes[0].set_ylabel('fraction'); axes[1].set_ylabel('effective mass (MeV)')
    axes[2].set_ylabel(r'$p_{Fn}-p_{Fp}-p_{F\ell}$ (MeV)'); axes[2].axhline(0, color='0.3', lw=.8)
    axes[2].set_xlabel(r'$n_B$ (fm$^{-3}$)'); axes[2].set_ylim(-150, 150)
    for a in axes:
        a.grid(alpha=.22); a.legend(fontsize=7.5, ncol=3)
    fig.suptitle('Cold npeμ matter used for the Urca comparison')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'eos_inputs_multi.{ext}', dpi=170)
    plt.close(fig)


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--points', type=int, default=120)
    a.add_argument('--eos', nargs='+', default=['DD2', 'BCP', 'BSk22'], choices=list(E.LOADERS))
    args = a.parse_args()
    results = {}
    for name in args.eos:
        rows = scan(name, args.points)
        results[name] = rows
        te, tm = thresholds(rows, 'gap_e'), thresholds(rows, 'gap_mu')
        st = stars(name)
        print(f"{name}: Mmax={st['Mmax']:.3f} nc(Mmax)={st['nc_Mmax']:.3f} nc(1.4)={st.get('nc_1.4', np.nan):.3f}; "
              f"M at DU_e onset {[round(mass_at(st, v), 3) for v in te]}, DU_mu {[round(mass_at(st, v), 3) for v in tm]}")
        print(f'{name}: {len(rows)} densities {rows[0]["nb"]:.3f}-{rows[-1]["nb"]:.3f}; '
              f'DU_e at {[round(v, 4) for v in te]}, DU_mu at {[round(v, 4) for v in tm]}; '
              f'max|beta res| {max(abs(r["beta_residual"]) for r in rows):.2e} MeV, '
              f'max|charge res| {max(abs(r["charge_residual"]) for r in rows):.1e}', flush=True)
    plot(results)
    print('outputs in', OUT)


if __name__ == '__main__':
    main()
