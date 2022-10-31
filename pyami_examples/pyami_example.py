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

def construct_example4():
    
    g = pyami.g_prod_t()

    alpha1 = pyami.VectorInt([1, 0, 0, 1, -1])
    alpha2 = pyami.VectorInt([0, 0, 0, 1, 0])
    alpha3 = pyami.VectorInt([1, 0, 0, 0, 0])
    alpha4 = pyami.VectorInt([1, 0, 0, 0, 0])
    alpha5 = pyami.VectorInt([0, 1, 0, 0, 0])
    alpha6 = pyami.VectorInt([-1, 1, 1, 0, 0])
    alpha7 = pyami.VectorInt([0, 0, 1, 0, 0])   

    epsilon1 = pyami.VectorInt([0,0,0,0,0,1,0])
    epsilon2 = pyami.VectorInt([0,0,0,0,1,0,0])
    epsilon3 = pyami.VectorInt([1,0,0,0,0,0,0])
    epsilon4 = pyami.VectorInt([1,0,0,0,0,0,0])
    epsilon5 = pyami.VectorInt([0,1,0,0,0,0,0])
    epsilon6 = pyami.VectorInt([0,0,0,1,0,0,0])
    epsilon7 = pyami.VectorInt([0,0,1,0,0,0,0])



    g1 = pyami.AmiBase.g_struct(epsilon1, alpha1)
    g2 = pyami.AmiBase.g_struct(epsilon2, alpha2)
    g3 = pyami.AmiBase.g_struct(epsilon3, alpha3)
    g4 = pyami.AmiBase.g_struct(epsilon4, alpha4)
    g5 = pyami.AmiBase.g_struct(epsilon5, alpha5)
    g6 = pyami.AmiBase.g_struct(epsilon6, alpha6)
    g7 = pyami.AmiBase.g_struct(epsilon7, alpha7)


    R0 = pyami.g_prod_t([g1, g2, g3, g4, g5, g6, g7])

    return R0

def construct_ext_example4():

    energy = pyami.VectorComplex([1, 1.1, 1.2, 1.31, 1.4, 0.01, 0.1])
    frequency = pyami.VectorComplex()
    for i in range(5):
        frequency.append(0+0j)
    
    frequency.append(0+math.pi*1j)
    beta = 1.0
    external = pyami.AmiBase.ami_vars(energy, frequency, beta)

    return external

def construct_example4():
    
    g = pyami.g_prod_t()

    alpha1 = pyami.VectorInt([1, 0, 0, 1, -1])
    alpha2 = pyami.VectorInt([0, 0, 0, 1, 0])
    alpha3 = pyami.VectorInt([1, 0, 0, 0, 0])
    alpha4 = pyami.VectorInt([1, 0, 0, 0, 0])
    alpha5 = pyami.VectorInt([0, 1, 0, 0, 0])
    alpha6 = pyami.VectorInt([-1, 1, 1, 0, 0])
    alpha7 = pyami.VectorInt([0, 0, 1, 0, 0])   

    epsilon1 = pyami.VectorInt([0,0,0,0,0,1,0])
    epsilon2 = pyami.VectorInt([0,0,0,0,1,0,0])
    epsilon3 = pyami.VectorInt([1,0,0,0,0,0,0])
    epsilon4 = pyami.VectorInt([1,0,0,0,0,0,0])
    epsilon5 = pyami.VectorInt([0,1,0,0,0,0,0])
    epsilon6 = pyami.VectorInt([0,0,0,1,0,0,0])
    epsilon7 = pyami.VectorInt([0,0,1,0,0,0,0])



    g1 = pyami.AmiBase.g_struct(epsilon1, alpha1)
    g2 = pyami.AmiBase.g_struct(epsilon2, alpha2)
    g3 = pyami.AmiBase.g_struct(epsilon3, alpha3)
    g4 = pyami.AmiBase.g_struct(epsilon4, alpha4)
    g5 = pyami.AmiBase.g_struct(epsilon5, alpha5)
    g6 = pyami.AmiBase.g_struct(epsilon6, alpha6)
    g7 = pyami.AmiBase.g_struct(epsilon7, alpha7)


    R0 = pyami.g_prod_t([g1, g2, g3, g4, g5, g6, g7])

    return R0

def construct_ext_example6():

    energy = pyami.VectorComplex([1,1.1,1.2,1.3,1.4,0, 0.1, 0.2, 0.3,0.4, 0.5])
    frequency = pyami.VectorComplex()
    for i in range(6):
        frequency.append(0+0j)
    
    frequency.append(0+math.pi*1j)
    beta = 1.0
    external = pyami.AmiBase.ami_vars(energy, frequency, beta)

    return external

def construct_example6():
    
    g = pyami.g_prod_t()

    alpha1=pyami.VectorInt([1,0,0,0,0,0,0])
    alpha2=pyami.VectorInt([0,1,0,0,0,0,0])
    alpha3=pyami.VectorInt([0,0,1,0,0,0,0])
    alpha4=pyami.VectorInt([0,0,0,1,0,0,0])
    alpha5=pyami.VectorInt([0,0,0,0,1,0,0])
    alpha6=pyami.VectorInt([0,0,0,0,0,1,0])
    alpha7=pyami.VectorInt([1,-1,0,0,0,1,0])
    alpha8=pyami.VectorInt([-1,1,0,0,0,0,1])
    alpha9=pyami.VectorInt([-1,1,0,0,0,0,1])
    alpha10=pyami.VectorInt([-1,1,-1,1,0,0,1])
    alpha11=pyami.VectorInt([1,0,0,0,-1,1,0])

    epsilon1=pyami.VectorInt([1,0,0,0,0,0,0,0,0,0,0])
    epsilon2=pyami.VectorInt([0,1,0,0,0,0,0,0,0,0,0])
    epsilon3=pyami.VectorInt([0,0,1,0,0,0,0,0,0,0,0])
    epsilon4=pyami.VectorInt([0,0,0,1,0,0,0,0,0,0,0])
    epsilon5=pyami.VectorInt([0,0,0,0,1,0,0,0,0,0,0])
    epsilon6=pyami.VectorInt([0,0,0,0,0,1,0,0,0,0,0])
    epsilon7=pyami.VectorInt([0,0,0,0,0,0,1,0,0,0,0])
    epsilon8=pyami.VectorInt([0,0,0,0,0,0,0,1,0,0,0])
    epsilon9=pyami.VectorInt([0,0,0,0,0,0,0,1,0,0,0])
    epsilon10=pyami.VectorInt([0,0,0,0,0,0,0,0,0,1,0])
    epsilon11=pyami.VectorInt([0,0,0,0,0,0,0,0,0,0,1])

    g1 = pyami.AmiBase.g_struct(epsilon1, alpha1)
    g2 = pyami.AmiBase.g_struct(epsilon2, alpha2)
    g3 = pyami.AmiBase.g_struct(epsilon3, alpha3)
    g4 = pyami.AmiBase.g_struct(epsilon4, alpha4)
    g5 = pyami.AmiBase.g_struct(epsilon5, alpha5)
    g6 = pyami.AmiBase.g_struct(epsilon6, alpha6)
    g7 = pyami.AmiBase.g_struct(epsilon7, alpha7)
    g8 = pyami.AmiBase.g_struct(epsilon8, alpha8)
    g9 = pyami.AmiBase.g_struct(epsilon9, alpha9)
    g10 = pyami.AmiBase.g_struct(epsilon10, alpha10)
    g11 = pyami.AmiBase.g_struct(epsilon11, alpha11)


    R0 = pyami.g_prod_t([g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11])

    return R0