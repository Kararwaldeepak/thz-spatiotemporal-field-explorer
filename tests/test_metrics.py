import unittest

from thz_field_analysis import calculate_scan_metrics


class TestScanMetrics(unittest.TestCase):
    def test_201_frames_10um_retroreflector(self):
        metrics = calculate_scan_metrics(
            stage_step_um=10.0,
            n_time_samples=201,
            delay_multiplier=2.0,
        )
        self.assertAlmostEqual(metrics.time_step_s * 1e15, 66.712819, places=5)
        self.assertAlmostEqual(metrics.frequency_resolution_hz / 1e12, 0.074575, places=5)
        self.assertAlmostEqual(metrics.nyquist_frequency_hz / 1e12, 7.494811, places=5)
        self.assertEqual(metrics.n_positive_frequency_frames, 101)
        self.assertAlmostEqual(metrics.stage_scan_span_um, 2000.0, places=8)

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            calculate_scan_metrics(0, 201, 2)
        with self.assertRaises(ValueError):
            calculate_scan_metrics(10, 1, 2)


if __name__ == "__main__":
    unittest.main()
