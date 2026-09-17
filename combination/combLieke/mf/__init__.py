"""TRExFitter MultiFit of Z -> ee, Z -> mumu and Z -> tau_h tau_h (CMS Open Data 2016 G+H).

Modules
    paths       repository and work-directory locations
    trexcfg     read, adapt and write the channels' TRExFitter configs; the MultiFit config
    ee_input    the z-ee histograms, taken from the TRExFitter job in z-ee/Zee_fit.tar.gz
    prediction  the aMC@NLO sigma(60 < m < 120 GeV) and its PDF / alpha_s / scale uncertainty
    results     parse the TRExFitter outputs into output/result.json
    plots       the CMS Open Data figures in output/plots/
"""
