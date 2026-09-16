#!/usr/bin/env python3
"""Prepare isolated TRExFitter configs, audit inputs, run h/w + MultiFit, export pb."""
import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import acceptance
import ee_archive
import histogram_audit

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
KINDS = {'Job', 'Fit', 'Region', 'Sample', 'NormFactor', 'Systematic'}


def read_config(path):
    """Read the simple HIST configs used here; reject unsupported block syntax."""
    blocks = []
    for number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(('%', '#')):
            continue
        # Preserve quoted ROOT labels containing # or %, remove unquoted comments.
        line = re.split(r'[%#](?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', line)[0].strip()
        key, sep, value = line.partition(':')
        if not sep:
            raise ValueError(f'{path}:{number}: unsupported syntax')
        value = value.strip().strip('"')
        if key in KINDS:
            blocks.append([key, value, {}])
        elif blocks:
            if key in blocks[-1][2]:
                raise ValueError(f'{path}:{number}: duplicate option {key}')
            blocks[-1][2][key] = value
        else:
            raise ValueError(f'{path}:{number}: option before block')
    return blocks


def render(blocks):
    # TRExFitter does not strip quotes for enums, numbers, booleans or lists.
    # Only human-readable labels/titles require quoting (notably ROOT's # markup).
    def value(key, val):
        return f'"{val}"' if ('Label' in key or 'Title' in key or key in {'Category', 'SubCategory'}) else str(val)
    return '\n\n'.join(f'{kind}: "{name}"\n' + '\n'.join(
        f'  {key}: {value(key, val)}' for key, val in opts.items())
        for kind, name, opts in blocks) + '\n'


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a finite positive number')
    return value


def prepare(manifest, *, run_dir=None, acceptance_factors=None, preliminary=False):
    run_dir = HERE if run_dir is None else Path(run_dir)
    acceptance_factors = acceptance_factors or dict.fromkeys(manifest["channels"], 1.0)
    if set(acceptance_factors) != {"mumu", "tautau"}:
        raise ValueError("Acceptance factors must specify mumu and tautau")
    for channel, factor in acceptance_factors.items():
        positive(factor, f"{channel}.acceptance_factor")
    acceptance_problems = acceptance.validate(manifest["acceptance"])
    reference = positive(manifest['reference_pb'], 'reference_pb')
    if set(manifest['channels']) != {'mumu', 'tautau'}:
        raise ValueError('Exactly mumu and tautau are required for this isolated run')
    generated = run_dir / 'generated'
    generated.mkdir(parents=True, exist_ok=True)
    problems, reviews, details, regions_seen = [], [], {}, set()
    multi = [['MultiFit', 'combLieke', dict(
        OutputDir=str(run_dir / 'results'), Combine='TRUE', Compare='FALSE',
        POIName='mu_comb', POIRange='0:3', DataName='obsData',
        Label='Z to ll: mumu + tautau', CmeLabel='13 TeV',
        LumiLabel='16.4 fb^{-1}', NumCPU='4')]]
    for channel, spec in manifest['channels'].items():
        source = REPO / spec['config']
        archive_info = None
        if spec.get('result_archive'):
            archive_info = ee_archive.read_archive(REPO / spec['result_archive'])
            (generated / f'{channel}_archive.json').write_text(json.dumps(archive_info, indent=2) + '\n')
        blocks = read_config(source)
        jobs = [b for b in blocks if b[0] == 'Job']
        if len(jobs) != 1 or jobs[0][2].get('ReadFrom') != 'HIST':
            raise ValueError(f'{channel}: expected one ReadFrom HIST job')
        job = jobs[0]
        oldpoi = job[2]['POI']
        norms = [b for b in blocks if b[0] == 'NormFactor' and b[1] == oldpoi]
        if len(norms) != 1:
            raise ValueError(f'{channel}: expected one POI NormFactor')
        signal = norms[0][2]['Samples']
        # Absolute paths make both individual and MultiFit reads independent of cwd.
        for kind, name, opts in blocks:
            if any(k in opts for k in ('Include', 'ReplacementFile', 'HistoPaths', 'HistoFiles')):
                raise ValueError(f'{channel}: includes/plural paths require an explicit adapter')
            if 'HistoPath' in opts:
                opts['HistoPath'] = str((source.parent / opts['HistoPath']).resolve())
            if kind == 'Systematic' and name in spec.get('nuisance_parameters', {}):
                opts['NuisanceParameter'] = spec['nuisance_parameters'][name]
            if kind == 'Systematic' and name == 'SigModel':
                opts['NuisanceParameter'] = f'SigModel_{channel}'
            if kind == 'Fit':
                opts['UseMinos'] = 'all' if opts.get('UseMinos', '').lower() == 'all' else 'mu_comb'
            if kind == 'Sample' and name in signal.split(','):
                opts['Type'] = 'SIGNAL'
        job[1] = f'comb_{channel}'
        job[2].update(POI='mu_comb', OutputDir=str(run_dir / 'results'),
                      HistoPath=str((REPO / spec['histogram_dir']).resolve()),
                      MCstatThreshold='0')
        # ee's bare SR is renamed; explicit HistoName remains h_mass.
        rename = {b[1]: f'{channel}_{b[1]}' for b in blocks
                  if b[0] == 'Region' and not b[1].startswith(channel + '_')}
        for b in blocks:
            if b[0] == 'Region':
                b[2].setdefault('HistoName', b[1])
                b[1] = rename.get(b[1], b[1])
                if b[1] in regions_seen:
                    raise ValueError(f'Duplicate region {b[1]}')
                regions_seen.add(b[1])
            if 'Regions' in b[2]:
                b[2]['Regions'] = ','.join(rename.get(x, x) for x in b[2]['Regions'].split(','))
        norms[0][1] = 'mu_comb'
        norms[0][2].update(Nominal='1', Min='0', Max='3', Title=f'sigma / {reference:g} pb')
        channel_ref = spec['sigma_reference_pb']
        if channel_ref is None:
            problems.append(f'{channel}: supply sigma_reference_pb for the agreed 60–120 GeV observable')
        else:
            positive(channel_ref, f'{channel}.sigma_reference_pb')
            factor = reference / channel_ref * acceptance_factors[channel]
            blocks.append(['NormFactor', f'reference_{channel}', dict(
                Nominal=str(factor), Min=str(factor), Max=str(factor),
                Constant='TRUE', Samples=signal)])
        if not spec['measurement_ready']:
            reviews.append(f'{channel}: measurement_ready=false; {spec["pending"]}')
        if spec['acceptance_in_likelihood']:
            problems.append(f'{channel}: remove in-likelihood acceptance variations before external refits to avoid double counting')
        if not any(b[0] == 'Systematic' for b in blocks):
            problems.append(f'{channel}: no systematic uncertainties in source config')
        files = set()
        for kind, name, opts in blocks:
            if kind == 'Sample':
                directory = Path(opts.get('HistoPath', job[2]['HistoPath']))
                filename = opts.get('HistoFile', job[2].get('HistoFile'))
                if not filename:
                    raise ValueError(f'{channel}/{name}: missing HistoFile')
                files.add(str(directory / (filename + '.root')))
        for filename in sorted(files):
            if not Path(filename).is_file():
                problems.append(f'{channel}: missing ROOT input {filename}')
        audit = histogram_audit.audit(blocks)
        problems.extend(f'{channel}: {problem}' for problem in audit['problems'])
        target = generated / f'{channel}.config'
        target.write_text('% Generated by combine.py; edit inputs.json or the source config.\n' + render(blocks))
        directory = run_dir / 'results' / job[1]
        multi.append(['Fit', job[1], dict(ConfigFile=str(target), Directory=str(directory), Label=channel)])
        details[channel] = dict(config=str(target), job=job[1], source=str(source),
                                root_files=sorted(files), sigma_reference_pb=channel_ref, archive_reference=archive_info,
                                histogram_audit=audit, systematics=[dict(name=b[1],
                                nuisance_parameter=b[2].get('NuisanceParameter', b[1]),
                                type=b[2]['Type'], samples=b[2].get('Samples', 'all'))
                                for b in blocks if b[0] == 'Systematic'])
    if not manifest['common_observable_and_overlap_validated']:
        reviews.append('Confirm common truth mass definition, lepton universality and event orthogonality')
    (generated / 'multifit.config').write_text(render(multi))
    blocking = problems + ([] if preliminary else reviews)
    status = dict(ready=not blocking, problems=blocking, technical_problems=problems,
                  review_pending=reviews, preliminary=preliminary, channels=details,
                  reference_pb=reference, observable=manifest['observable'],
                  run_dir=str(run_dir), acceptance_factors=acceptance_factors,
                  acceptance_ready=not acceptance_problems, acceptance_problems=acceptance_problems)
    (run_dir / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
    return status


def ensure_workspace_alias(status, info):
    """Require the fitted-region workspace; allBinsFitRegions restores excluded bins."""
    job = info['job']
    directory = Path(status['run_dir']) / 'results' / job / 'RooStats'
    conventional = directory / f'{job}_combined_{job}_model.root'
    if not conventional.is_file():
        raise ValueError(f'Fitted-region workspace missing after hw: {conventional}; do not substitute an all-bin workspace')
    return conventional


def run_fit(status):
    if not status['ready']:
        raise ValueError('Inputs incomplete: see status.json. No fit was started.')
    exe = shutil.which('trex-fitter')
    if not exe:
        raise ValueError('trex-fitter missing: source the repository setup.sh (v1.8.0 required)')
    run_dir = Path(status['run_dir'])
    logs = run_dir / 'logs'
    logs.mkdir(exist_ok=True)
    def run(action, config, label):
        print(f'trex-fitter {action} {config}', flush=True)
        with (logs / f'{label}.log').open('w') as output:
            subprocess.run([exe, action, str(config)], cwd=run_dir,
                           stdout=output, stderr=subprocess.STDOUT, check=True)
    # Always rebuild: never combine stale workspaces with an old normalization/NP map.
    for channel, info in status['channels'].items():
        run('hw', info['config'], channel)
        ensure_workspace_alias(status, info)
    run('mwf', run_dir / 'generated/multifit.config', 'multifit')
    return export_result(status)


def export_result(status):
    path = Path(status['run_dir']) / 'results/combLieke/Fits/combLieke.txt'
    match = re.search(r'^mu_comb\s+([-+\d.eE]+)\s+\+([-+\d.eE]+)\s+-([-+\d.eE]+)',
                      path.read_text(), re.MULTILINE)
    if not match:
        raise ValueError(f'No mu_comb result in {path}')
    mu, up, down = map(float, match.groups())
    if not all(math.isfinite(x) for x in (mu, up, down)) or min(up, down) < 0:
        raise ValueError('Invalid POI result')
    ref = status['reference_pb']
    result = dict(observable=status['observable'], reference_pb=ref,
                  mu_comb=mu, sigma_pb=mu * ref, error_up_pb=up * ref,
                  error_down_pb=down * ref, fit_file=str(path),
                  validation='Inspect logs, fit status, covariance quality, pulls and POI boundaries before quoting.')
    result['preliminary'] = status.get('preliminary', False)
    result['review_pending'] = status.get('review_pending', [])
    result['uncertainty_scope'] = 'in-fit only; acceptance fixed'
    result['acceptance_factors'] = status['acceptance_factors']
    (Path(status['run_dir']) / 'fit_result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def run_combination(manifest, status, *, nominal_only=False):
    deferred = manifest['acceptance']['method'] == 'postfit_deferred'
    if not nominal_only and not deferred and not status['acceptance_ready']:
        raise ValueError('External acceptance inputs incomplete: ' + '; '.join(status['acceptance_problems']))
    # A separate file distinguishes a complete combination from an in-fit-only run.
    summary_path = HERE / 'result.json'
    if summary_path.exists():
        summary_path.unlink()
    nominal = run_fit(status)
    if nominal_only or deferred:
        nominal['acceptance_status'] = 'pending_postfit'
        nominal['acceptance_note'] = ('Nominal channel acceptance is included in the reference normalization. '
            'Acceptance uncertainty remains separate and must be propagated after the MultiFit, '
            'following channel prescriptions; no combined acceptance error or total has been assigned.')
        (Path(status['run_dir']) / 'fit_result.json').write_text(json.dumps(nominal, indent=2) + '\n')
        return nominal
    variations = {}
    for source, direction, factors in acceptance.scenarios(manifest['acceptance']):
        run_dir = HERE / 'variations' / source['name'] / direction
        shifted = prepare(manifest, run_dir=run_dir, acceptance_factors=factors)
        variations.setdefault(source['name'], {})[direction] = run_fit(shifted)
    result = acceptance.summarize(nominal, variations, manifest['acceptance'])
    summary_path.write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preliminary', action='store_true', help='allow pending scientific reviews; retain all technical input checks and label results preliminary')
    parser.add_argument('--inputs', type=Path, default=HERE / 'inputs.json')
    parser.add_argument('--run', action='store_true', help='run nominal MultiFit; acceptance follows the manifest policy (default: deferred post-fit)')
    parser.add_argument('--nominal-only', action='store_true', help='with --run, fit at fixed acceptance without external errors')
    args = parser.parse_args()
    if args.nominal_only and not args.run:
        parser.error('--nominal-only requires --run')
    try:
        manifest = json.loads(args.inputs.read_text())
        status = prepare(manifest, preliminary=args.preliminary)
        print('Prepared mumu and tautau channel configs and generated/multifit.config')
        if args.preliminary:
            for review in status['review_pending']:
                print('PRELIMINARY REVIEW PENDING:', review)
        for problem in status['problems']:
            print('NOMINAL INPUT NEEDED:', problem)
        for problem in status['acceptance_problems']:
            print('POST-FIT ACCEPTANCE PENDING:' if manifest['acceptance']['method'] == 'postfit_deferred'
                  else 'EXTERNAL ACCEPTANCE INPUT NEEDED:', problem)
        if args.run:
            result = run_combination(manifest, status, nominal_only=args.nominal_only)
            print(json.dumps(result, indent=2))
        else:
            print('Nominal fit ready' if status['ready'] else 'Nominal fit requires the inputs above.')
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(2, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
