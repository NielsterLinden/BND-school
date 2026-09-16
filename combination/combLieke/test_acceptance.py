import copy
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import acceptance
import combine


def model():
    return dict(method='external_refits', inputs_complete=True,
                combine_with_infit_quadrature=True, sources=[dict(
                    name='Example', description='Synthetic fully correlated response for tests only',
                    relative_shifts={ch: dict(up=0.1, down=-0.1) for ch in acceptance.CHANNELS})])


class AcceptanceTests(unittest.TestCase):
    def test_validation(self):
        for value in (-1, float('nan'), float('inf'), True, None):
            m = model()
            m['sources'][0]['relative_shifts']['ee']['up'] = value
            with self.assertRaises(ValueError):
                acceptance.validate(m)
        m = model()
        del m['sources'][0]['relative_shifts']['ee']
        with self.assertRaises(ValueError):
            acceptance.validate(m)
        m = model()
        m['sources'].append(copy.deepcopy(m['sources'][0]))
        with self.assertRaises(ValueError):
            acceptance.validate(m)

    def test_incomplete_acceptance_does_not_block_nominal(self):
        manifest = json.loads((combine.HERE / 'inputs.json').read_text())
        with tempfile.TemporaryDirectory(dir=combine.HERE) as tmp:
            with patch.object(combine, 'HERE', Path(tmp)):
                status = combine.prepare(manifest)
                self.assertFalse(status['acceptance_ready'])
                self.assertFalse(any('in-likelihood acceptance variations' in p or 'Acceptance inputs are incomplete' in p for p in status['problems']))
                with patch.object(combine, 'run_fit', return_value={'sigma_pb': 2000}) as run:
                    result = combine.run_combination(manifest, status)
                    self.assertEqual(result['sigma_pb'], 2000)
                    self.assertEqual(result['acceptance_status'], 'pending_postfit')
                    self.assertNotIn('approximate_total', result)
                    self.assertFalse((Path(tmp) / 'result.json').exists())
                    run.assert_called_once()
                manifest['acceptance']['method'] = 'external_refits'
                with patch.object(combine, 'run_fit') as run:
                    with self.assertRaises(ValueError):
                        combine.run_combination(manifest, status)
                    run.assert_not_called()

    def test_refit_pipeline_inverse_scaling_and_isolation(self):
        manifest = json.loads((combine.HERE / 'inputs.json').read_text())
        manifest['acceptance'] = model()
        manifest['channels']['ee']['sigma_reference_pb'] = 1900
        with tempfile.TemporaryDirectory(dir=combine.HERE) as tmp:
            with patch.object(combine, 'HERE', Path(tmp)):
                status = combine.prepare(manifest)
                nominal_text = Path(status['channels']['mumu']['config']).read_text()
                seen = []
                def fake_fit(state):
                    # Controlled fixed-yield likelihood: sigma * A is constant.
                    seen.append(Path(state['run_dir']))
                    for ch, spec in state['channels'].items():
                        blocks = combine.read_config(Path(spec['config']))
                        nf = next(b[2] for b in blocks if b[1] == f'reference_{ch}')
                        expected = state['reference_pb'] / spec['sigma_reference_pb'] * state['acceptance_factors'][ch]
                        self.assertAlmostEqual(float(nf['Nominal']), expected)
                    return dict(sigma_pb=2000 / state['acceptance_factors']['mumu'], error_up_pb=20, error_down_pb=18)
                with patch.object(combine, 'run_fit', side_effect=fake_fit):
                    result = combine.run_combination(manifest, status)
                self.assertEqual(len(set(seen)), 3)
                self.assertEqual(Path(status['channels']['mumu']['config']).read_text(), nominal_text)
                self.assertAlmostEqual(result['acceptance']['error_up_pb'], 2000 / .9 - 2000)
                self.assertAlmostEqual(result['acceptance']['error_down_pb'], 2000 - 2000 / 1.1)
                self.assertAlmostEqual(result['approximate_total']['error_up_pb'], math.hypot(20, 2000 / .9 - 2000))
                self.assertTrue((Path(tmp) / 'result.json').is_file())

    def test_independent_modes_envelopes_and_no_total(self):
        m = model()
        m['combine_with_infit_quadrature'] = False
        m['sources'].append(dict(name='Other', description='Independent synthetic mode',
                                  relative_shifts={ch: dict(up=0, down=0) for ch in acceptance.CHANNELS}))
        result = acceptance.summarize(dict(sigma_pb=100, error_up_pb=2, error_down_pb=2), {
            'Example': {'up': {'sigma_pb': 103}, 'down': {'sigma_pb': 99}},
            'Other': {'up': {'sigma_pb': 104}, 'down': {'sigma_pb': 102}}}, m)
        self.assertEqual(result['acceptance']['error_up_pb'], 5)
        self.assertEqual(result['acceptance']['error_down_pb'], 1)
        self.assertNotIn('approximate_total', result)

    def test_channel_anticorrelation_and_zero_response(self):
        m = model()
        m['sources'][0]['relative_shifts']['ee'] = dict(up=-.1, down=.1)
        m['sources'][0]['relative_shifts']['tautau'] = dict(up=0, down=0)
        cases = list(acceptance.scenarios(m))
        self.assertEqual(cases[0][2], dict(mumu=1.1, ee=.9, tautau=1))
        self.assertEqual(cases[1][2], dict(mumu=.9, ee=1.1, tautau=1))


if __name__ == '__main__':
    unittest.main()
