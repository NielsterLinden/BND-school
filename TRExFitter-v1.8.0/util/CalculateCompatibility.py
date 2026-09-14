"""Calculate compatibility between different hypotheses, see also the TRExFitter website
https://trexfitter-docs.web.cern.ch/trexfitter-docs/latest/model_studies/compatibility/"""
from __future__ import print_function

import sys

if sys.argv[-1] == __file__ or len(sys.argv) != 4:
    print("run like this:")
    print("python", __file__, "1mu_NLL 2mu_NLL #dof")
    print("example: python", __file__, "11.8 11.3 1")
    sys.exit()


null_NLL = float(sys.argv[1])  # null hypothesis NLL
alt_NLL = float(sys.argv[2])  # alternative hypothesis NLL, smaller than null_NLL
ndof = int(sys.argv[3])

try:
    # use ROOT if available
    from ROOT import TMath

    p_val = TMath.Prob(2 * (null_NLL - alt_NLL), ndof)
    print("The compatibility is " + str(p_val * 100) + "%.")

except ImportError:
    # if no ROOT is available, try scipy
    from scipy.stats import chi2

    p_val = chi2.sf(2 * (null_NLL - alt_NLL), ndof)  # sf is the same as 1-cdf
    print("The compatibility is " + str(p_val * 100) + "%.")
