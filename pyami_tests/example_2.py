import pyami

ami = pyami.AmiBase()

# external variables
energy = pyami.VectorDouble([-4, 0.1, -1])
freq = [3.14j, 3.14j, 3.14j]
beta = 1.00
avars = pyami.AmiBase.ami_vars(energy, freq, beta)

# second order self energy setup
a1 = pyami.VectorInt([1, 0, 0])
a2 = pyami.VectorInt([0, 1, 0])
a3 = pyami.VectorInt([-1, 1, 1])

e1 = pyami.VectorInt([1, 0, 0])
e2 = pyami.VectorInt([0, 1, 0])
e3 = pyami.VectorInt([0, 0, 1])

g1 = pyami.AmiBase.g_struct(e1, a1)
g2 = pyami.AmiBase.g_struct(e2, a2)
g3 = pyami.AmiBase.g_struct(e3, a3)

R0 = [g1, g2, g3]

# These are opaque types so that they are mutable in c++ construct function.
S_array = pyami.S_t()
P_array = pyami.P_t()
R_array = pyami.R_t()

E_REG = 0 # numberical regulator
N_INT = 2 # number of matsubara sums
test_amiparms = pyami.AmiBase.ami_parms(N_INT, E_REG)

# now construct!
ami.construct(test_amiparms, R0, R_array, P_array, S_array)

# Get integrand for these specific external params
ans = ami.evaluate(test_amiparms, R_array, P_array, S_array, avars)

print(ans)
