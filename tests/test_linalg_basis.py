import unittest
from celltools.linalg.basis import Basis, Vector, LinearAlgebraError
import numpy as np
from numpy.testing import assert_array_equal

class TestBasis(unittest.TestCase):
    def test_valid_initialization(self):
        v1, v2, v3 = np.array([1,0,1]), np.array([0,1,0]), np.array([1,2,3])
        base = Basis(v1,v2,v3, offset=(5,6,7))
        assert_array_equal(base.offset, np.array([5,6,7]), 
                        "Basis initialize: Offset not set correctly")
        assert_array_equal(base.basis, np.array([v1,v2,v3]),
                        "Basis initialize: Basis not set correctly")

    def test_linear_independence(self):
        with self.assertRaises(LinearAlgebraError) as context:
            Basis([1,0,0],[0,1,0],[1,1,0])
        self.assertEqual(str(context.exception), "vectors are not linear independent", "Basis linear independence: linear dependence not caught" )

class TestVector(unittest.TestCase):

    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)
        self.v1, self.v2, self.v3, self.v4 = np.array([1,0,1]), np.array([0,1,0]), np.array([1,2,3]), np.array([2,4,6])
        self.base1 = Basis(self.v1,self.v2,self.v3,offset=(5,6,7))
        self.base2 = Basis(self.v1,self.v2,self.v4)

    def test_valid_initialization(self):
        vec = Vector([1,0,0], self.base1)
        assert_array_equal(vec.global_coord, self.v1, "Vector initialize: Vector not initialized correctly")
    
    def test_vector_sum(self):
        vec1 = Vector([1,0,0], self.base1)
        vec2 = Vector([0,1,0], self.base1)
        vec3 = Vector([-1,2,3], self.base2)
        assert_array_equal((vec1 + vec2).global_coord, np.array([1,1,1]),"Vector sum: Sum result incorrect")
        assert_array_equal((vec1 - vec2).global_coord, np.array([1,-1,1]),"Vector sum: Subtraction result incorrect")
        

        with self.assertRaises(LinearAlgebraError) as context:
            vec1 + vec3
        self.assertEqual(str(context.exception),"basis do not match" ,"Vector sum: Different basis in sum not caught")
        
        with self.assertRaises(LinearAlgebraError) as context:
            vec1 - vec3
        self.assertEqual(str(context.exception),"basis do not match" ,"Vector sum: Different basis in subtraction not caught")


    def test_scalar_multiplication(self):
        vec1 = Vector([1,0,0], self.base1)

        assert_array_equal((3.5*vec1).global_coord, 3.5*self.v1, "Vector scalar multiplication: Multiplication result incorrect")

    def test_parallel(self):
        vec1 = Vector([0,0,1], self.base1)
        vec2 = Vector([0,0,1], self.base2)
        vec3 = Vector([1,1,1], self.base1)
        self.assertEqual(vec1.parallel(vec2), True, "Vector parallel: Parallel vectors not detected")
        self.assertEqual(vec1.parallel(vec3), False, "Vector parallel: Not parallel vectors not detected")

if __name__ == "__main__":
    unittest.main()

