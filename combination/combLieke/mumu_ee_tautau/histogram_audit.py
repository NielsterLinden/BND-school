"""Read-only validation of the single-file HIST conventions used by these channels."""
from contextlib import ExitStack
from pathlib import Path


def applies(options, key, name):
    values = {x.strip() for x in options.get(key, 'all').split(',')}
    excluded = {x.strip() for x in options.get('Exclude' + key, '').split(',')}
    return ('all' in values or name in values) and name not in excluded


def audit(blocks):
    report = dict(checked=False, nominal_histograms=0, variation_histograms=0,
                  problems=[], warnings=[])
    try:
        import numpy as np
        import uproot
    except ImportError:
        report['problems'].append('Histogram contents not checked: source setup.sh for uproot and rerun preparation')
        return report
    job = next(b[2] for b in blocks if b[0] == 'Job')
    regions = [b for b in blocks if b[0] == 'Region']
    samples = [b for b in blocks if b[0] == 'Sample']
    systematics = [b for b in blocks if b[0] == 'Systematic' and b[2]['Type'] == 'HISTO']
    with ExitStack() as stack:
        files = {}
        def read(path, name, is_mc, edges=None):
            try:
                if path not in files:
                    files[path] = stack.enter_context(uproot.open(path))
                h = files[path][name]
                if not h.classname.startswith('TH1'):
                    raise ValueError('expected one-dimensional TH1')
                values, actual_edges = h.to_numpy()
                variances = h.variances(flow=True)
                if not np.all(np.isfinite(h.values(flow=True))):
                    raise ValueError('nonfinite contents')
                if variances is None or not np.all(np.isfinite(variances)) or np.any(variances < 0):
                    raise ValueError('invalid variances')
                if is_mc and len(h.member('fSumw2')) != len(values) + 2:
                    raise ValueError('missing explicit MC Sumw2')
                if edges is not None and not np.array_equal(edges, actual_edges):
                    raise ValueError('binning differs from nominal/region')
                if np.any(values < 0):
                    report['warnings'].append(f'{path.name}:{name}: negative bins; review with channel DropBins/rebinning')
                return actual_edges
            except (OSError, KeyError, ValueError) as error:
                report['problems'].append(f'{path.name}:{name}: {str(error).splitlines()[0]}')
                return None
        for _, region, ro in regions:
            region_edges = None
            for _, sample, so in samples:
                if not applies(so, 'Regions', region):
                    continue
                path = Path(so.get('HistoPath', job['HistoPath'])) / (so.get('HistoFile', job.get('HistoFile', '')) + '.root')
                nominal = so.get('HistoName', ro.get('HistoName', region)) + so.get('HistoNameSuff', '')
                edges = read(path, nominal, so['Type'] != 'DATA', region_edges)
                report['nominal_histograms'] += 1
                if edges is not None:
                    region_edges = edges
                if so['Type'] == 'DATA':
                    continue
                for _, syst, opts in systematics:
                    if not applies(opts, 'Samples', sample) or not applies(opts, 'Regions', region):
                        continue
                    for direction in ('Up', 'Down'):
                        # Fail closed if a future delivery needs more than this adapter supports.
                        if any(k in opts for k in ('HistoPath' + direction, 'HistoFile' + direction)):
                            report['problems'].append(f'{syst}: variation file/path override requires explicit audit support')
                            continue
                        name = opts.get('HistoName' + direction, nominal + opts.get('HistoNameSuf' + direction, ''))
                        read(path, name, True, edges)
                        report['variation_histograms'] += 1
    report['checked'] = True
    return report
