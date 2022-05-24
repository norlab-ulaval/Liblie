import unittest
import scipy as sc
import numpy as np

import liblie as ll

class TestStringMethods(unittest.TestCase):
    def setUp(self):
        self.M2 = np.array([[1,2],[3,4]])
        self.M3 = np.array([[1,2,3],[4,5,6],[7,8,9]])
        
        self.coords_so3 = np.array([.1,.2,.3])
        self.so3_M3 = np.array([[0,-.3,.2],[.3,0,-.1],[-.2,.1,0]]) #norm of coordinate <pi

        self.coords_se3 = np.array([1,2,3,.1,.2,.3])
        self.se3_M4 = np.zeros((4,4))
        self.se3_M4[:3,:3] = self.so3_M3
        self.se3_M4[:3,3] = np.array([1,2,3])
        self.SO3_M3 = sc.linalg.expm(self.so3_M3)
        self.SE3_M4 = sc.linalg.expm(self.se3_M4)

        self.se3_M4_curly = np.zeros((6,6))
        self.se3_M4_curly[:3,:3] = self.so3_M3
        self.se3_M4_curly[3:,3:] = self.so3_M3
        self.se3_M4_curly[:3,3:] = np.array([[0,-3,2],[3,0,-1],[-2,1,0]])

        self.a1 = np.array([.2,.2,0])
        self.a2 = np.array([0,0,.01])
        self.C1 = sc.linalg.expm(ll.wedge_so3(self.a1))
        self.C2 = sc.linalg.expm(ll.wedge_so3(self.a2))
        self.C = self.C1 @ self.C2

        self.b1 = np.array([1,2,3,.2,.2,0])
        self.b2 = np.array([.01,0,0,0,0,.01])
        self.T1 = sc.linalg.expm(ll.wedge_se3(self.b1))
        self.T2 = sc.linalg.expm(ll.wedge_se3(self.b2))
        self.T = self.T1 @ self.T2

    def assertArrayEqual(self, x, y,rtol=1e-7, atol=1e-7): #wrapper
        #np.testing.assert_array_equal(x,y)
        np.testing.assert_allclose(x,y, rtol=rtol, atol=atol)

#------------------------------------------------------------------
    def test_expm(self):
        self.assertArrayEqual(ll.expm(self.M2), sc.linalg.expm(self.M2))
        self.assertArrayEqual(ll.expm(self.M3), sc.linalg.expm(self.M3))

        self.assertArrayEqual(ll.expm_so3(self.so3_M3), sc.linalg.expm(self.so3_M3))
        self.assertArrayEqual(ll.expm_se3(self.se3_M4), sc.linalg.expm(self.se3_M4))

        self.assertArrayEqual(ll.expm_so3(np.zeros((3,3))), np.eye(3,3))
        self.assertArrayEqual(ll.expm_se3(np.zeros((4,4))), np.eye(4,4))

    def test_logm(self):
        self.assertArrayEqual(ll.logm(self.M2), sc.linalg.logm(self.M2))
        self.assertArrayEqual(ll.logm(self.M3), sc.linalg.logm(self.M3))

        self.assertArrayEqual(ll.logm_so3(self.SO3_M3), self.so3_M3)
        self.assertArrayEqual(ll.logm_se3(self.SE3_M4), self.se3_M4)

        self.assertArrayEqual(ll.logm_so3(np.eye(3,3)), np.zeros((3,3)))
        self.assertArrayEqual(ll.logm_se3(np.eye(4,4)), np.zeros((4,4)))

    def test_wedge_vee(self):
        self.assertArrayEqual(ll.wedge_so3(self.coords_so3), self.so3_M3)
        self.assertArrayEqual(self.coords_so3, ll.vee_so3(self.so3_M3))

        self.assertArrayEqual(ll.wedge_se3(self.coords_se3), self.se3_M4)
        self.assertArrayEqual(self.coords_se3, ll.vee_se3(self.se3_M4))
        self.assertArrayEqual(ll.curly_wedge_se3(self.coords_se3), self.se3_M4_curly)

    def test_jacobians(self):
        #exp(a1^)exp(a2^) ~= [a1 + Jl(a1)^-1 * a2]^ for small a2
        self.assertArrayEqual(
                ll.vee_so3(sc.linalg.logm(self.C)),
                self.a1 + ll.jac_inv_so3(self.a1)@self.a2, rtol=1e-3, atol=1e-2)
        self.assertArrayEqual(
                ll.jac_inv_so3(self.a1)@ll.jac_so3(self.a1),
                np.eye(3,3))

        #se3 jacobians are only approximations, so we need to raise the tolerance
        self.assertArrayEqual(
                ll.vee_se3(sc.linalg.logm(self.T)),
                self.b1 + ll.jac_inv_se3(self.b1)@self.b2, rtol=1e-3, atol=1e-1)
        self.assertArrayEqual(
                ll.jac_inv_se3(self.b1)@ll.jac_se3(self.b1),
                np.eye(6,6), rtol=1e-2, atol=5e-1)

        



if __name__ == '__main__':
    unittest.main()
