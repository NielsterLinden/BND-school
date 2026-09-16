#!/usr/bin/env python3
"""Publication-style comparisons using the saved nominal MultiFit result."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / 'combination'))
from comb import plots
import matplotlib.pyplot as plt

result = json.loads((HERE / 'fit_result.json').read_text())
manifest = json.loads((HERE / 'inputs.json').read_text())
out = HERE / 'output/plots'
pred = sum(c['sigma_reference_pb'] for c in manifest['channels'].values()) / 2
combined = dict(kind='this', label=r'combined MultiFit: $\mu\mu + \tau\tau$',
                value=result['sigma_pb'], down=result['error_down_pb'], up=result['error_up_pb'],
                colour=plots.C_COMB)

def reference(label, value, error, citation):
    return dict(kind='pub', label=label, value=value, down=error, up=error,
                colour=plots.C_GREY, citation=citation)

references = [reference(*r) for r in plots.REFERENCES]
mumu = json.loads((REPO / 'z-mumu/fit/results/zmumu_fit_result.json').read_text())
tau = json.loads((REPO / 'z-tautau/fit/results/ztautau_fit_result.json').read_text())
channel_rows = []
for name, value, down, up in [
    ('mumu', mumu['mu']*manifest['channels']['mumu']['sigma_reference_pb'],
     mumu['mu_err_down']*manifest['channels']['mumu']['sigma_reference_pb'],
     mumu['mu_err_up']*manifest['channels']['mumu']['sigma_reference_pb']),
    ('tautau', tau['poi_value']*manifest['channels']['tautau']['sigma_reference_pb'],
     tau['poi_err_down']*manifest['channels']['tautau']['sigma_reference_pb'],
     tau['poi_err_up']*manifest['channels']['tautau']['sigma_reference_pb'])]:
    channel_rows.append(dict(kind='this', label=plots._label(name) + ' standalone',
                             value=value, down=down, up=up, colour=plots.COLOUR[name]))
    channel_rows += [reference(label, v, e, cite) for n,label,v,e,cite in plots.CHANNEL_REFERENCES if n==name]
channel_rows += [combined] + references

def draw(rows, filename):
    rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(8.6, 1.65 + 0.55*len(rows)))
    ax.axvline(pred, color=plots.C_PRED, lw=1.3, ls='--', label=f'aMC@NLO reference {pred:.0f} pb')
    for i, row in enumerate(rows):
        own = row['kind']=='this'
        ax.errorbar(row['value'], i, xerr=[[row['down']], [row['up']]],
                    fmt='o' if own else 's', color=row['colour'], ms=7 if own else 5.5,
                    capsize=4, lw=2.4 if own else 1.6)
        text = (fr"{row['value']:.1f} $^{{+{row['up']:.1f}}}_{{-{row['down']:.1f}}}$ pb" if own
                else fr"{row['value']:.0f} $\pm$ {row['up']:.0f} pb")
        ax.text(row['value'], i+0.23, text, ha='center', fontsize=8.5,
                color=row['colour'], fontweight='bold' if own else 'normal')
    ax.set_yticks(range(len(rows)), [('this work: ' if r['kind']=='this' else '')+r['label'] for r in rows],fontsize=8.5)
    for tick, row in zip(ax.get_yticklabels(), rows):
        if row['kind']=='this': tick.set_fontweight('bold')
    ax.set_ylim(-0.7,len(rows)-0.25)
    lo=min(r['value']-r['down'] for r in rows); hi=max(r['value']+r['up'] for r in rows)
    pad=0.12*(hi-lo);ax.set_xlim(min(lo-pad,pred-35),max(hi+pad,pred+35))
    ax.set_xlabel(r'$\sigma(pp\to Z/\gamma^*\to\ell\ell)$, $60<m<120$ GeV [pb]')
    ax.set_title('CMS Open Data 2016, 16.4 fb$^{-1}$ — preliminary\n'
                 'This work: full in-fit uncertainty; acceptance uncertainty excluded',loc='left',fontsize=10)
    ax.grid(axis='y',visible=False); ax.legend(fontsize=9,loc='lower right')
    fig.text(0.01, -0.01, 'Published points: total uncertainties; ATLAS 66–116 GeV scaled by 1.01425.\n'
             'Reference values and conversion inherited from combination/comb/plots.py.', fontsize=8,color=plots.C_GREY)
    plots._save(fig,out/filename)

draw(channel_rows,'comparison_channels')
draw([combined]+references,'comparison')
(out/'comparison_inputs.json').write_text(json.dumps(dict(combined=combined,channel_comparison_rows=channel_rows,
    reference_pb=pred,uncertainty_scope='this work: in-fit only; published: total',
    reference_source='combination/comb/plots.py'),indent=2)+'\n')
print(out)
