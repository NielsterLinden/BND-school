"""Z -> mu+ mu- fiducial cross-section measurement on CMS 2016 Open Data.

The package is deliberately small and explicit: every physics choice lives in
`config.py`, every reusable routine in one of the four modules below.

    config   selection cuts, luminosity, paths, binning   (edit this first)
    io       file discovery, chunked reading, golden-JSON lumi mask
    objects  muon selection, FSR recovery, dimuon kinematics
    hists    histogram booking, saving and CMS-style plotting
    stats    efficiency intervals and uncertainty propagation

See docs/ for the physics documentation and README.md for how to run.
"""

__version__ = "1.0.0"
