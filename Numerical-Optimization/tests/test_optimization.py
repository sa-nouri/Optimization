import unittest
import numpy as np
from src.levenberg_marquardt import LevenbergMarquardt
from src.simulated_annealing import SimulatedAnnealing
from src.conj_grad_linear import ConjugateGradient

class TestOptimizationAlgorithms(unittest.TestCase):
    def setUp(self):
        # Test function: f(x) = x^2 + 2x + 1
        self.test_function = lambda x: x**2 + 2*x + 1
        self.test_gradient = lambda x: 2*x + 2
        self.test_hessian = lambda x: 2

    def test_levenberg_marquardt(self):
        lm = LevenbergMarquardt(self.test_function, self.test_gradient, self.test_hessian)
        x0 = np.array([1.0])
        result = lm.optimize(x0)
        self.assertAlmostEqual(result[0], -1.0, places=5)  # Minimum is at x = -1

    def test_simulated_annealing(self):
        sa = SimulatedAnnealing(self.test_function)
        x0 = np.array([1.0])
        result = sa.optimize(x0)
        self.assertAlmostEqual(result[0], -1.0, places=2)  # Less precise due to stochastic nature

    def test_conjugate_gradient(self):
        # Test with a simple quadratic function
        A = np.array([[2, 0], [0, 2]])
        b = np.array([-2, -2])
        cg = ConjugateGradient()
        x0 = np.array([1.0, 1.0])
        result = cg.solve(A, b, x0)
        expected = np.array([-1.0, -1.0])
        np.testing.assert_array_almost_equal(result, expected)

if __name__ == '__main__':
    unittest.main() 
