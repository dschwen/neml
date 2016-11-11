import sys
sys.path.append('..')

from neml import hardening
import unittest

from common import *

import numpy as np
import numpy.linalg as la
import numpy.random as ra


class CommonHardening(object):
  """
    Tests that can apply to all hardening rules
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_gradient(self):
    dfn = lambda x: self.model.q(x, self.T)
    ngrad = differentiate(dfn, self.hist_trial)
    grad = self.model.dq_da(self.hist_trial, self.T)
    self.assertTrue(np.allclose(ngrad, grad))

class TestLinearIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.abs(ra.random((1,)))
    self.T = 300.0

    self.model = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.s0, self.s0))
    self.assertTrue(np.isclose(self.model.K, self.K))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      np.array([-self.s0 - self.K * self.hist_trial[0]])))

class TestVoceIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.R = 100.0
    self.d = 10.0

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.abs(ra.random((1,)))
    self.T = 300.0

    self.model = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.s0, self.s0))
    self.assertTrue(np.isclose(self.model.R, self.R))
    self.assertTrue(np.isclose(self.model.d, self.d))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      np.array([-self.s0 - self.R * (1 - np.exp(-self.d*self.hist_trial[0]))])))

class TestLinearKinematicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.H = 1000.0

    self.hist0 = np.zeros((6,))
    
    self.hist_trial = ra.random((6,)) * 100
    self.hist_trial = self.hist_trial - np.array([1,1,1,0,0,0]) * sum(self.hist_trial[:3]) / 3.0
    self.T = 300.0

    self.model = hardening.LinearKinematicHardeningRule(self.H)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.H, self.H))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      -self.hist_trial * self.H))

class TestCombinedHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0
    self.H = 1000.0

    self.hist0 = np.zeros((7,))
    self.hist_trial = ra.random((7,)) * 100
    self.hist_trial[1:] = make_dev(self.hist_trial[1:])
    self.T = 300.0

    self.iso = hardening.LinearIsotropicHardeningRule(self.s0, self.K)
    self.kin = hardening.LinearKinematicHardeningRule(self.H)

    self.model = hardening.CombinedHardeningRule(self.iso, self.kin)

  def test_relation(self):
    sb = np.zeros((7,))
    sb[0] = self.iso.q(self.hist_trial[0:1], self.T)
    sb[1:] = self.kin.q(self.hist_trial[1:], self.T)

    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), sb))

class CommonNonAssociative(object):
  """
    Common tests for non-associative hardening rules
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertEqual(self.model.ninter, self.conform)
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def gen_stress(self):
    s = ra.random((6,))
    s = (1.0 - 2.0 * s) * 125.0
    return s

  def test_dq(self):
    a = self.gen_hist()

    dq_model = self.model.dq_da(a, self.T)
    dfn = lambda x: self.model.q(x, self.T)
    dq_num = differentiate(dfn, a)

    self.assertTrue(np.allclose(dq_model, dq_num))

  def test_dh_ds(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_ds(s, a, self.T)
    dfn = lambda x: self.model.h(x, a, self.T)
    dh_num = differentiate(dfn, s)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3))

  def test_dh_da(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_da(s, a, self.T)
    dfn = lambda x: self.model.h(s, x, self.T)
    dh_num = differentiate(dfn, a)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3))

class TestChaboche(unittest.TestCase, CommonNonAssociative):
  """
    Chaboche model with arbitrary kinematic hardening
  """
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.n = 4
    self.cs = ra.random((self.n,)) * 10.0
    self.rs = ra.random((self.n,)) * 10.0

    self.iso = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

    self.model = hardening.Chaboche(self.iso, self.cs, self.rs)

    self.hist0 = np.zeros((1 + self.n*6,))
    self.conform = 7

    self.T = 300.0

  def gen_hist(self):
    hist = ra.random((1 + self.n*6,))
    hist[1:] = (1.0 - 2.0 * hist[1:]) * 100.0
    for i in range(self.n):
      hist[1+i*6:1+(i+1)*6] = make_dev(hist[1+i*6:1+(i+1)*6])
    return hist

  def test_properties(self):
    self.assertEqual(self.n, self.model.n)
    self.assertTrue(np.allclose(self.model.c, self.cs))
    self.assertTrue(np.allclose(self.model.r, self.rs))

  def test_q(self):
    h = self.gen_hist()
    q_model = self.model.q(h, self.T)
    q_exact = np.zeros((7,))
    q_exact[0] = self.iso.q(h[0:1], self.T)
    q_exact[1:] = sum(h[1+i*6:1+(i+1)*6] for i in range(self.n))
    self.assertTrue(np.allclose(q_model, q_exact))

  def test_h(self):
    alpha = self.gen_hist()
    s = self.gen_stress()
    sdev = make_dev(s)
    X = sum(alpha[1+i*6:1+(i+1)*6] for i in range(self.n))
    n = (sdev+X) / la.norm(sdev+X)

    h_model = self.model.h(s, alpha, self.T)

    h_exact = np.zeros((self.model.nhist,))
    h_exact[0] = np.sqrt(2.0/3.0)
    for i in range(self.n):
      h_exact[1+i*6:1+(i+1)*6] = -self.cs[i] * self.rs[i] * (n + alpha[1+i*6:1+(i+1)*6] / self.rs[i])
    
    self.assertTrue(np.allclose(h_model, h_exact))