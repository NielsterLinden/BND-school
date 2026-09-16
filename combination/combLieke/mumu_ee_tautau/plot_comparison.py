#!/usr/bin/env python3
"""Plot saved three-channel MultiFit and channel results against verified paper values."""
import json
import math
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO/'combination'))
from comb import plots as style
read = lambda p: json.loads(p.read_text())
result = read(HERE/'fit_result.json')
manifest = read(HERE/'inputs.json')
papers = read(HERE/'paper_references.json')
previous = read(HERE.parent/'mumu_tautau/fit_result.json')
out = HERE/'output/plots'
out.mkdir(parents=True, exist_ok=True)

def own(label, value, down, up, colour, source):
    assert all(math.isfinite(x) for x in [value,down,up]) and min(down,up)>0
    return dict(kind='this',label=label,value_pb=value,error_down_pb=down,error_up_pb=up,
                colour=colour,source=source,uncertainty_scope='in-fit only; acceptance uncertainty excluded')

def from_result(label, r, colour, source):
    return own(label,r['sigma_pb'],r['error_down_pb'],r['error_up_pb'],colour,source)

m = read(REPO/'z-mumu/fit/results/zmumu_fit_result.json')
t = read(REPO/'z-tautau/fit/results/ztautau_fit_result.json')
e = read(HERE/'ee_fit_result.json')
refm = manifest['channels']['mumu']['sigma_reference_pb']
reft = manifest['channels']['tautau']['sigma_reference_pb']
channels = {
    'mumu': own(r'This work: $Z\to\mu\mu$',m['mu']*refm,m['mu_err_down']*refm,m['mu_err_up']*refm,style.C_MUMU,'z-mumu/fit/results/zmumu_fit_result.json'),
    'ee': from_result(r'This work: $Z\to ee$ (new inputs)',e,style.C_EE,'ee_fit_result.json'),
    'tautau': own(r'This work: $Z\to\tau_h\tau_h$ (Tight)',t['poi_value']*reft,t['poi_err_down']*reft,t['poi_err_up']*reft,style.C_TAUTAU,'z-tautau/fit/results/ztautau_fit_result.json')}
combined = from_result(r'This work: $\mu\mu+ee+\tau\tau$ MultiFit',result,style.C_COMB,'fit_result.json')
prev = from_result(r'Previous: $\mu\mu+\tau\tau$ MultiFit',previous,'#577D29','../mumu_tautau/fit_result.json')

def published(r):
    channel = {'mumu':r'$Z\to\mu\mu$', 'ee':r'$Z\to ee$', 'tautau':r'$Z\to\tau\tau$', 'combined':r'$ee+\mu\mu$'}[r['channel']]
    label = f"{r['experiment']} {r['year']}: {channel}"
    if r['experiment']=='ATLAS': label += ' *'
    return dict(kind='pub',label=label,value_pb=r['value_pb'],error_down_pb=r['error_pb'],error_up_pb=r['error_pb'],colour=style.C_GREY,source=r['source'],uncertainty_scope='published total',reference=r)

rows=[]
for channel in ['mumu','ee','tautau']:
    rows.append(channels[channel])
    rows.extend(published(r) for r in papers['references'] if r['channel']==channel)
rows += [prev,combined]
rows.extend(published(r) for r in papers['references'] if r['channel']=='combined')
pred=sum(c['sigma_reference_pb'] for c in manifest['channels'].values())/3

def draw(rows,filename):
    fig,ax=plt.subplots(figsize=(9.2,2.2+len(rows)*0.47))
    fig.subplots_adjust(left=0.36,right=0.78,bottom=0.21,top=0.85)
    ax.axvline(pred,color=style.C_PRED,lw=1.3,ls='--',label=f'aMC@NLO reference {pred:.0f} pb')
    for i,r in enumerate(rows):
        is_this=r['kind']=='this'
        if r['label']==combined['label']:ax.axhspan(i-0.43,i+0.43,color='#edf1f1',zorder=0)
        ax.errorbar(r['value_pb'],i,xerr=[[r['error_down_pb']],[r['error_up_pb']]],fmt='o' if is_this else 's',
                    color=r['colour'],ms=6.5 if is_this else 5.2,capsize=4,lw=2.2 if is_this else 1.6)
        text=(fr"{r['value_pb']:.1f} $^{{+{r['error_up_pb']:.1f}}}_{{-{r['error_down_pb']:.1f}}}$" if is_this
              else fr"{r['value_pb']:.1f} $\pm$ {r['error_up_pb']:.1f}")
        ax.text(1.035,i,text,transform=ax.get_yaxis_transform(),va='center',fontsize=9,color=r['colour'],fontweight='bold' if is_this else 'normal')
    ax.set_yticks(range(len(rows)),[r['label'] for r in rows],fontsize=9)
    for tick,r in zip(ax.get_yticklabels(),rows):
        if r['kind']=='this':tick.set_fontweight('bold')
    ax.set_ylim(len(rows)-0.4,-0.65)
    lo=min(r['value_pb']-r['error_down_pb'] for r in rows);hi=max(r['value_pb']+r['error_up_pb'] for r in rows)
    pad=(hi-lo)*.12;ax.set_xlim(min(lo-pad,pred-35),max(hi+pad,pred+35))
    ax.set_xlabel(r'$\sigma(pp\to Z/\gamma^*\to\ell\ell)$ [pb]',fontsize=11)
    ax.grid(axis='y',visible=False)
    ax.text(1.035,1.025,'Cross section [pb]',transform=ax.transAxes,fontsize=9)
    fig.text(.03,.965,'Z cross sections at 13 TeV',fontsize=16,fontweight='bold')
    fig.text(.03,.929,'CMS Open Data 2016, 16.4 fb$^{-1}$ | 60 < m < 120 GeV | Preliminary',fontsize=11)
    fig.text(.03,.895,'This work: in-fit uncertainties; acceptance uncertainty excluded. Published points: total uncertainties.',fontsize=9)
    fig.text(.03,.075,'* ATLAS: 66–116 GeV values × 1.01425 (analysis mass-window conversion).\n'
             'Papers: ATLAS 1603.09222; CMS 2408.03744 (published 2025) and 1801.03535.\n'
             'Standalone and combined results share data; these rows are not independent measurements.',fontsize=8.5,color='#555555')
    ax.legend(loc='lower right',fontsize=8)
    style._save(fig,out/filename)

draw(rows,'comparison_channels')
draw([prev,combined]+[published(r) for r in papers['references'] if r['channel']=='combined'],'comparison')
(out/'comparison_inputs.json').write_text(json.dumps(dict(rows=rows,reference_pb=pred,fit_gof_pvalue=result.get('reported_gof_pvalue'),paper_references='../../paper_references.json'),indent=2)+'\n')
print(out)
