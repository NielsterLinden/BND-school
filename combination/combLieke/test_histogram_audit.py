import copy
import unittest
from pathlib import Path

import combine
import histogram_audit


class HistogramAuditTests(unittest.TestCase):
    def setUp(self):
        try:
            import uproot
        except ImportError:
            self.skipTest('source setup.sh for uproot')
        self.blocks = combine.read_config(combine.REPO / 'z-ee/fit.config')
        self.blocks[0][2]['HistoPath'] = str(combine.REPO / 'datasets/z-ee')

    def test_delivery_and_missing_variation(self):
        report = histogram_audit.audit(self.blocks)
        self.assertTrue(report['checked'])
        self.assertEqual(report['problems'], [])
        self.assertEqual(report['nominal_histograms'], 6)
        self.assertEqual(report['variation_histograms'], 44)
        self.assertEqual(report['warnings'], [])
        changed = copy.deepcopy(self.blocks)
        next(b[2] for b in changed if b[0] == 'Systematic' and b[1] == 'PDF')['HistoNameUp'] = 'missing_PDFUp'
        report = histogram_audit.audit(changed)
        self.assertTrue(any('missing_PDFUp' in p for p in report['problems']))

    def test_sample_and_region_filtering(self):
        self.assertTrue(histogram_audit.applies({}, 'Samples', 'signal'))
        self.assertFalse(histogram_audit.applies({'Samples': 'background'}, 'Samples', 'signal'))
        self.assertFalse(histogram_audit.applies({'ExcludeRegions': 'SR'}, 'Regions', 'SR'))
