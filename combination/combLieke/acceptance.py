"""External acceptance variations: independent modes, correlated responses within a mode."""
import math
import re

CHANNELS = {'mumu', 'ee', 'tautau'}


def validate(model):
    if model.get('method') == 'postfit_deferred':
        return ['Acceptance uncertainty deferred until the nominal MultiFit succeeds; follow channel post-fit prescriptions']
    if model.get('method') != 'external_refits':
        raise ValueError('acceptance.method must be postfit_deferred or external_refits')
    if not isinstance(model.get('combine_with_infit_quadrature'), bool):
        raise ValueError('acceptance.combine_with_infit_quadrature must be boolean')
    sources = model.get('sources')
    if not isinstance(sources, list):
        raise ValueError('acceptance.sources must be a list')
    seen = set()
    for source in sources:
        name = source['name']
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', name) or name in seen:
            raise ValueError('Acceptance source names must be unique safe identifiers')
        seen.add(name)
        if not isinstance(source.get('description'), str) or not source['description'].strip():
            raise ValueError(f'{name}: document the source and correlation assumption')
        responses = source['relative_shifts']
        if set(responses) != CHANNELS:
            raise ValueError(f'{name}: specify mumu, ee and tautau explicitly (zero if unaffected)')
        for channel, response in responses.items():
            if set(response) != {'up', 'down'}:
                raise ValueError(f'{name}/{channel}: supply up and down signed relative A shifts')
            for shift in response.values():
                if isinstance(shift, bool) or not isinstance(shift, (int, float)) or not math.isfinite(shift) or shift <= -1:
                    raise ValueError(f'{name}/{channel}: shifts must be finite numbers greater than -1')
    problems = []
    if model.get('inputs_complete') is not True:
        problems.append('Acceptance inputs are incomplete: supply signed responses and correlation definitions')
    if not sources:
        problems.append('No external acceptance sources supplied')
    return problems


def scenarios(model):
    problems = validate(model)
    if problems:
        raise ValueError('; '.join(problems))
    for source in model['sources']:
        for direction in ('up', 'down'):
            yield source, direction, {
                channel: 1 + response[direction]
                for channel, response in source['relative_shifts'].items()}


def summarize(nominal, variations, model):
    """Envelope each independent mode, then quadrature; retain signed fit responses."""
    validate(model)
    nominal_sigma = nominal['sigma_pb']
    components = {}
    for source in model['sources']:
        name = source['name']
        shifts = {direction: variations[name][direction]['sigma_pb'] - nominal_sigma
                  for direction in ('up', 'down')}
        components[name] = dict(
            description=source['description'], relative_shifts=source['relative_shifts'],
            shift_up_variation_pb=shifts['up'], shift_down_variation_pb=shifts['down'],
            error_up_pb=max(0, *shifts.values()),
            error_down_pb=max(0, *(-x for x in shifts.values())))
    acc = {f'error_{direction}_pb': math.sqrt(sum(
        component[f'error_{direction}_pb'] ** 2 for component in components.values()))
        for direction in ('up', 'down')}
    result = dict(nominal=nominal, sigma_pb=nominal_sigma,
                  in_fit={key: nominal[key] for key in ('error_up_pb', 'error_down_pb')},
                  acceptance=dict(**acc, sources=components), variation_results=variations,
                  assumptions=[
                      'Different acceptance sources are independent modes; channels within each mode vary together.',
                      'Each shifted fit profiles the original in-fit nuisances; A–C theory correlations are not propagated jointly.',
                      'Source errors use the envelope of the two signed shifts, combined in quadrature across modes.'])
    if model['combine_with_infit_quadrature']:
        result['approximate_total'] = {
            f'error_{direction}_pb': math.hypot(nominal[f'error_{direction}_pb'], acc[f'error_{direction}_pb'])
            for direction in ('up', 'down')}
        result['assumptions'].append('The approximate total treats external acceptance and in-fit errors as independent.')
    return result
