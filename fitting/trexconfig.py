"""Render a TRExFitter (v1.8.0) configuration file from python structures.

The channel code builds a list of blocks and calls `render`:

    from fitting import trexconfig as tc
    blocks = [tc.job("zmumu", ...), tc.fit(...), tc.region("mumu_SR", ...), ...]
    Path("fit/zmumu.config").write_text(tc.render(blocks))

A block is (kind, name, ordered list of (key, value)). Values: str, int, float, bool,
list/tuple (comma-joined), or None (skipped). Strings are quoted unless they look like a
keyword/number/list (TRExFitter accepts both; quoting protects labels with spaces or '#').
Option names are those of TRExFitter v1.8.0 (`jobSchema.config` in its source tree).
"""

from __future__ import annotations

_KEYWORDS = {"TRUE", "FALSE", "HIST", "NTUP", "SPLUSB", "BONLY", "CRSR", "CRONLY", "SIGNAL",
             "BACKGROUND", "DATA", "GHOST", "CONTROL", "VALIDATION", "HISTO", "OVERALL",
             "SHAPE", "STAT", "ONESIDED", "TWOSIDED", "ABSMEAN", "MAXIMUM", "NONE", "ALL",
             "MERGE", "SYSTS", "GAMMAS", "ASIMOV", "STANDARD", "POISSON", "GAUSSIAN", "NOCRASH"}


def _fmt(value):
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (int, float)):
        return repr(value) if isinstance(value, float) else str(value)
    if isinstance(value, (list, tuple)):
        return ",".join(_fmt(v).strip('"') for v in value)
    s = str(value)
    if s in _KEYWORDS or s.replace(".", "", 1).replace("-", "", 1).isdigit():
        return s
    if "," in s and " " not in s and "#" not in s:       # sample/region lists
        return s
    return f'"{s}"'


def block(kind: str, name: str, **opts):
    """Generic block; keyword order is preserved."""
    return (kind, name, [(k, v) for k, v in opts.items() if v is not None])


def job(name, **opts):
    return block("Job", name, **opts)


def fit(name="fit", **opts):
    return block("Fit", name, **opts)


def region(name, **opts):
    return block("Region", name, **opts)


def sample(name, **opts):
    return block("Sample", name, **opts)


def normfactor(name, **opts):
    return block("NormFactor", name, **opts)


def overall_syst(name, up, down, samples, category, title=None, regions=None, subcategory=None):
    """Normalisation-only systematic (relative, e.g. up=0.012, down=-0.012)."""
    return block("Systematic", name, Title=title or name, Type="OVERALL", OverallUp=float(up),
                 OverallDown=float(down), Samples=samples, Regions=regions, Category=category,
                 SubCategory=subcategory)


def histo_syst(name, samples, category, title=None, regions=None, symmetrisation="TWOSIDED",
               smoothing=None, one_sided=False, drop_norm=None, subcategory=None):
    """Histogram systematic read as <region>__<sample>__<name>Up/Down."""
    opts = dict(Title=title or name, Type="HISTO", HistoNameSufUp=f"__{name}Up")
    if not one_sided:
        opts["HistoNameSufDown"] = f"__{name}Down"
    opts.update(Samples=samples, Regions=regions, Symmetrisation="ONESIDED" if one_sided else symmetrisation,
                Smoothing=smoothing, DropNorm=drop_norm, Category=category, SubCategory=subcategory)
    return block("Systematic", name, **opts)


def render(blocks, header: str = "") -> str:
    lines = []
    if header:
        lines += [f"% {line}" for line in header.rstrip().splitlines()]
        lines.append("")
    for kind, name, opts in blocks:
        lines.append(f'{kind}: "{name}"')
        for key, value in opts:
            lines.append(f"  {key}: {_fmt(value)}")
        lines.append("")
    return "\n".join(lines) + "\n"
