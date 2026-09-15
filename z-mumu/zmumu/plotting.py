"""Data / MC comparison plots from output/v2/histograms.pkl (stack + ratio, MC-stat band)."""

from __future__ import annotations

import numpy as np

from . import config, hists, histograms as H

STACK_ORDER = ["Fakes", "WJets", "WW", "WZ", "ZZ", "SingleTop", "TTbar", "DYee", "DYtautau", "DYmumu"]
COLOURS = {"DYmumu": "#f2b134", "DYtautau": "#a04cb0", "DYee": "#4c9ee0", "TTbar": "#d62728", "SingleTop": "#e07b7b",
           "WW": "#2ca6a4", "WZ": "#4dc0be", "ZZ": "#7fd4d2", "Fakes": "#9a9a9a", "WJets": "#8c6d4f"}
LABELS = {"DYmumu": r"Z/$\gamma^*\to\mu\mu$", "DYtautau": r"Z/$\gamma^*\to\tau\tau$", "DYee": r"Z/$\gamma^*\to ee$",
          "TTbar": r"t$\bar{t}$", "SingleTop": "tW", "WW": "WW", "WZ": "WZ", "ZZ": "ZZ", "Fakes": "non-prompt (FF)",
          "WJets": r"W+jets (jet $\to$ e)"}
XLABELS = {"mass_fit": r"$m_{\mu\mu}$ [GeV]", "mass_fine": r"$m_{\mu\mu}$ [GeV]", "pt1": r"leading muon $p_T$ [GeV]",
           "pt2": r"subleading muon $p_T$ [GeV]", "eta1": r"leading muon $\eta$", "eta2": r"subleading muon $\eta$",
           "phi1": r"leading muon $\phi$", "zpt": r"$p_T^{\mu\mu}$ [GeV]", "zy": r"$y^{\mu\mu}$", "npv": "good primary vertices",
           "met": r"PF $p_T^{miss}$ [GeV]", "puppimet": r"Puppi $p_T^{miss}$ [GeV]", "njet": "lepton-cleaned jets ($p_T$ > 30 GeV)", "iso1": "leading muon rel. iso.",
           "iso2": "subleading muon rel. iso.", "nfsr": "FSR photons", "pt_el": r"electron $p_T$ [GeV]"}


def stack_plot(hall, region, var, filename, fakes_fine=None, logy=False, title=""):
    import matplotlib.pyplot as plt
    import mplhep as hep
    e = H.edges(region, var)
    c = 0.5 * (e[1:] + e[:-1])
    data = hall.get(f"Data|{region}|{var}|nominal")
    if data is None:
        return None
    comps, labels, cols = [], [], []
    tot = np.zeros(len(c)); tot_var = np.zeros(len(c))
    for smp in STACK_ORDER:
        k = f"{smp}|{region}|{var}|nominal"
        v = hall.get(k)
        if smp == "Fakes":
            if region == "SR" and fakes_fine is not None and var in ("mass_fit", "mass_fine"):
                f = np.asarray(fakes_fine)
                v = np.clip(f.reshape(-1, len(f) // len(c)).sum(axis=1), 0, None)
            else:
                continue
            var_v = np.zeros_like(v)
        elif v is None or np.sum(v) <= 0:
            continue
        else:
            var_v = hall[k + "|w2"]
        comps.append(v); labels.append(LABELS[smp]); cols.append(COLOURS[smp])
        tot += v; tot_var += var_v
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(9, 9.5), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})
    # linear axes: plot in thousands when the scale is large (keeps the offset text off the CMS label)
    unit = 1.0
    if not logy and max(tot.max(), data.max()) >= 2e4:
        unit = 1e3
    hep.histplot([c_ / unit for c_ in comps], bins=e, ax=ax, stack=True, histtype="fill", label=labels, color=cols, edgecolor="black", linewidth=0.4)
    ax.fill_between(e, np.append(tot - np.sqrt(tot_var), 0) / unit, np.append(tot + np.sqrt(tot_var), 0) / unit, step="post",
                    facecolor="none", hatch="////", edgecolor="grey", linewidth=0, label="MC stat.")
    hep.histplot(data / unit, bins=e, ax=ax, histtype="errorbar", color="black", label="Data", markersize=4, yerr=np.sqrt(np.maximum(data, 0)) / unit)
    gev = " GeV" if var in ("mass_fit", "mass_fine", "pt1", "pt2", "zpt", "met", "pt_el") else ""
    ax.set_ylabel(f"Events / {e[1]-e[0]:g}{gev}" + (r" [$\times 10^3$]" if unit != 1.0 else ""))
    tot, data = tot / unit, data / unit
    if logy:
        ax.set_yscale("log"); ax.set_ylim(max(tot.min() * 0.1, 0.5), tot.max() * 30)
    else:
        ax.set_ylim(0, max(tot.max(), data.max()) * 1.35)
    ax.legend(fontsize=11, ncol=2, loc="upper right")
    hists._decorate(ax, lumi_fb=config.LUMI_PB_NORMTAG / 1000)
    hists._title(ax, title or {"SR": "signal region", "SS": "same-sign region", "CRemu": r"e$\mu$ region (OS), all MC",
                               "SSemu": r"e$\mu$ region (SS), all MC"}[region])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(tot > 0, data / tot, np.nan)
        err = np.where(tot > 0, np.sqrt(np.maximum(data * unit, 0)) / unit / tot, np.nan)
        band = np.where(tot > 0, np.sqrt(tot_var) / unit / tot, 0)
    rax.fill_between(e, np.append(1 - band, 1), np.append(1 + band, 1), step="post", facecolor="none", hatch="////", edgecolor="grey", linewidth=0)
    rax.errorbar(c, ratio, yerr=err, fmt="o", color="black", markersize=3)
    rax.axhline(1.0, linestyle="--", color="grey")
    xlabel = XLABELS.get(var, var)
    if region in H.EMU_REGIONS and var in ("mass_fit", "mass_fine"):
        xlabel = r"$m_{e\mu}$ [GeV]"
    rax.set_ylim(0.8, 1.2); rax.set_ylabel("Data / pred."); rax.set_xlabel(xlabel)
    return hists.save_fig(fig, filename)


def all_plots(hall, fakes_res=None):
    fine = np.array(fakes_res["templates"]["nominal"]) if fakes_res else None
    made = []
    for region, vars_ in H.REGION_VARS.items():
        for var in vars_:
            for logy in ((False, True) if var in ("mass_fit", "mass_fine", "pt1", "pt2", "zpt", "met", "puppimet") else (False,)):
                name = f"datamc_{region}_{var}{'_log' if logy else ''}.png"
                if stack_plot(hall, region, var, name, fine, logy=logy):
                    made.append(name)
    return made
