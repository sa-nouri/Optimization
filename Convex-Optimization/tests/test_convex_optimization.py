import unittest
import numpy as np
import cvxpy as cp

class TestConvexOptimization(unittest.TestCase):
    def setUp(self):
        # Simple linear programming problem
        self.x = cp.Variable(2)
        self.objective = cp.Minimize(self.x[0] + self.x[1])
        self.constraints = [
            self.x[0] >= 0,
            self.x[1] >= 0,
            self.x[0] + self.x[1] >= 1
        ]
        self.prob = cp.Problem(self.objective, self.constraints)

    def test_linear_programming(self):
        # Test basic linear programming
        self.prob.solve()
        self.assertAlmostEqual(self.prob.value, 1.0)
        self.assertAlmostEqual(self.x[0].value, 0.5)
        self.assertAlmostEqual(self.x[1].value, 0.5)

    def test_quadratic_programming(self):
        # Test quadratic programming
        x = cp.Variable(2)
        objective = cp.Minimize(cp.sum_squares(x))
        constraints = [x[0] + x[1] == 1]
        prob = cp.Problem(objective, constraints)
        
        prob.solve()
        self.assertAlmostEqual(prob.value, 0.5)
        self.assertAlmostEqual(x[0].value, 0.5)
        self.assertAlmostEqual(x[1].value, 0.5)

    def test_semidefinite_programming(self):
        # Test semidefinite programming
        X = cp.Variable((2, 2), symmetric=True)
        objective = cp.Minimize(cp.trace(X))
        constraints = [
            X >> 0,  # X is positive semidefinite
            X[0, 0] + X[1, 1] == 1
        ]
        prob = cp.Problem(objective, constraints)
        
        prob.solve()
        self.assertAlmostEqual(prob.value, 1.0)
        # The optimal X should be a 2x2 identity matrix scaled by 0.5
        self.assertAlmostEqual(X.value[0, 0], 0.5)
        self.assertAlmostEqual(X.value[1, 1], 0.5)
        self.assertAlmostEqual(X.value[0, 1], 0.0)
        self.assertAlmostEqual(X.value[1, 0], 0.0)

if __name__ == '__main__':
    unittest.main() 
