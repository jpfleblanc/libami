import pyami
import math
import time

def construct_example2():

    g = pyami.g_prod_t()

    alpha1 = pyami.VectorInt([1, 0, 0])
    alpha2 = pyami.VectorInt([0, 1, 0])
    alpha3 = pyami.VectorInt([-1, 1, 1])

    epsilon1 = pyami.VectorInt([1, 0, 0])
    epsilon2 = pyami.VectorInt([0, 1, 0])
    epsilon3 = pyami.VectorInt([0, 0, 1])

    g1 = pyami.AmiBase.g_struct(epsilon1, alpha1)
    g2 = pyami.AmiBase.g_struct(epsilon2, alpha2)
    g3 = pyami.AmiBase.g_struct(epsilon3, alpha3)

    R0 = pyami.g_prod_t([g1, g2, g3])

    return R0

def construct_ext_example2():

    energy = pyami.VectorComplex([-4, 0.1, -1])
    frequency = pyami.VectorComplex()
    for i in range(2):
        frequency.append(0+0j)
    
    frequency.append(0+math.pi*1j)
    beta = 1.0
    external = pyami.AmiBase.ami_vars(energy, frequency, beta)

    return external