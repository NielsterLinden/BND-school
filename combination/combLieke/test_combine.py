import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import combine


class CombinationTests(unittest.TestCase):
    def test_parser_preserves_root_labels(self):
        with tempfile.TemporaryDirectory(dir=combine.HERE) as folder:
            path = Path(folder) / 'test.config'
            path.write_text('Job: "test"\n Label: "#mu#mu" % comment\n ReadFrom: HIST\n')
            self.assertEqual(combine.read_config(path)[0][2]['Label'], '#mu#mu')

    def test_generation_and_normalization_closure(self):
        manifest = json.loads((combine.HERE / 'inputs.json').read_text())
        manifest['channels']['ee']['sigma_reference_pb'] = 1990.0
        with tempfile.TemporaryDirectory(dir=combine.HERE) as folder:
            with patch.object(combine, 'HERE', Path(folder)):
                status = combine.prepare(manifest)
                self.assertFalse(status['ready'])
                configs = Path(folder) / 'generated'
                multi = (configs / 'multifit.config').read_text()
                self.assertEqual(multi.count('\nFit:'), 3)
                for name, ch in manifest['channels'].items():
                    blocks = combine.read_config(configs / f'{name}.config')
                    factor = next(b[2] for b in blocks if b[1] == f'reference_{name}')
                    self.assertEqual(factor['Constant'], 'TRUE')
                    scale = float(factor['Nominal'])
                    self.assertAlmostEqual(scale * ch['sigma_reference_pb'] / manifest['reference_pb'], 1)
                    for b in blocks:
                        if b[0] == 'Systematic' and b[1] == 'SigModel':
                            self.assertEqual(b[2]['NuisanceParameter'], f'SigModel_{name}')
                            self.assertEqual(b[2]['HistoNameSufUp'], '__SigModelUp')
                ee_blocks = combine.read_config(configs / 'ee.config')
                lumi = next(b for b in ee_blocks if b[0] == 'Systematic' and b[1] == 'LUMI')
                self.assertEqual(lumi[2]['NuisanceParameter'], 'Lumi')
                self.assertEqual(len([b for b in ee_blocks if b[0] == 'Systematic']), 11)
                self.assertFalse(any('missing ROOT input' in x for x in status['problems']))
                with patch.object(combine.subprocess, 'run') as run:
                    with self.assertRaises(ValueError):
                        combine.run_fit(status)
                    run.assert_not_called()

    def test_dropbins_workspace_never_replaced_by_all_bins(self):
        with tempfile.TemporaryDirectory(dir=combine.HERE) as folder:
            directory = Path(folder) / 'results/comb_tautau/RooStats'
            directory.mkdir(parents=True)
            status, info = {'run_dir': folder}, {'job': 'comb_tautau'}
            all_bins = directory / 'comb_tautau_allBinsFitRegions_combined_comb_tautau_model.root'
            all_bins.write_bytes(b'all bins, including excluded bins')
            with self.assertRaises(ValueError):
                combine.ensure_workspace_alias(status, info)
            fitted = directory / 'comb_tautau_combined_comb_tautau_model.root'
            fitted.write_bytes(b'only fitted bins')
            self.assertEqual(combine.ensure_workspace_alias(status, info).read_bytes(), b'only fitted bins')

    def test_ee_archive_systematic_model(self):
        report = combine.ee_archive.read_archive(combine.REPO / 'z-ee/Zee_fit.tar.gz')
        self.assertAlmostEqual(report['poi_value'], 1.05087)
        self.assertIn('HistoSys', report['systematic_nodes'])
        self.assertIn('OverallSys', report['systematic_nodes'])
        self.assertGreater(report['error_decomposition']['SYST_ERROR']['symmetric'], 0)
        self.assertEqual(report['luminosity_model']['LumiRelErr'], '0')
        self.assertEqual(len(report['samples']), 5)
        self.assertTrue(all(s['mc_stat'] for s in report['samples']))
        self.assertTrue(report['workspace_available'])
        self.assertIsNone(report['sigma_reference_60_120_pb'])

    def test_preliminary_keeps_technical_checks(self):
        manifest = json.loads((combine.HERE / 'inputs.json').read_text())
        with tempfile.TemporaryDirectory(dir=combine.HERE) as folder:
            with patch.object(combine.histogram_audit, 'audit', return_value={'problems': []}):
                status = combine.prepare(manifest, run_dir=folder, preliminary=True)
                self.assertTrue(status['ready'])
                self.assertTrue(status['review_pending'])
                manifest['channels']['ee']['sigma_reference_pb'] = None
                status = combine.prepare(manifest, run_dir=folder, preliminary=True)
                self.assertFalse(status['ready'])
                self.assertTrue(any('sigma_reference_pb' in p for p in status['problems']))

    def test_render_trex_typed_options(self):
        text = combine.render([['Job', 'test', {'ReadFrom': 'HIST', 'NumCPU': '4',
            'Combine': 'TRUE', 'Samples': 'a,b', 'Label': '#mu#mu'}]])
        self.assertIn('ReadFrom: HIST', text)
        self.assertIn('NumCPU: 4', text)
        self.assertIn('Samples: a,b', text)
        self.assertIn('Label: "#mu#mu"', text)

    def test_invalid_reference(self):
        for value in (0, -1, float('nan'), float('inf'), None, True):
            with self.assertRaises(ValueError):
                combine.positive(value, 'reference')

    def test_missing_channel_rejected(self):
        manifest = json.loads((combine.HERE / 'inputs.json').read_text())
        del manifest['channels']['ee']
        with self.assertRaises(ValueError):
            combine.prepare(manifest)


if __name__ == '__main__':
    unittest.main()
