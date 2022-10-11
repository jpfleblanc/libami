import pyami

ami = pyami.AmiBase()

energy = pyami.VectorDouble([-4, 0.1, -1])
freq = [3.14j, 3.14j, 3.14j]
beta = 1.0

avars = pyami.AmiBase.ami_vars(energy, freq, beta)

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


# might run into trouble here, vectors in c++ are intialized to be push_backed on to
# This does not play nice in python :(
S_array = pyami.S_t()#[[pyami.VectorDouble()]]#[[pyami.VectorDouble()]] # inner most should be doubles
P_array = pyami.P_t()#[[[]]]#[[[[pyami.AmiBase.pole_struct()]]]] # inner most should be pole_structs
R_array = pyami.R_t()#[[[]]]#[[[pyami.AmiBase.g_struct()]]] # inner most should be g_structs

E_REG = 0
N_INT = 2
test_amiparms = pyami.AmiBase.ami_parms(N_INT, E_REG)

# now construct!
print("Before the construction: S_t, P_t and R_t are:")
print(S_array)
print(P_array)
print(R_array)
ami.construct(test_amiparms, R0, R_array, P_array, S_array)
print("After the construction: S_t, P_t and R_t are:")
print(S_array)
print(P_array)
print(R_array)

ans = ami.evaluate(test_amiparms, R_array, P_array, S_array, avars)

print(ans)
# segfault 11 :(



