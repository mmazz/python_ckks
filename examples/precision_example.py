import sys
from ckks import scheme
import time
import numpy as np
import mpmath as mp
arithmetic64 = False
# q0=scale=P=1024, level=2
"""
La multiplicacion con nivelacion de level funciona con
scale = 2**9
level = 4
q0 = 2**4
P = scale**level*q0
"""

if len(sys.argv) > 1:
    logM =  int(sys.argv[1])
    logp =  int(sys.argv[2])
    logq =  int(sys.argv[3])
else:
    logM =  3
    logp =  9
    logq =  4
sparce_ternary = False
M = 2**logM
roundCS = 2
h = 2
scale = 2**logp # el p del paper.
level = 3
q0 = 2**logq
P = scale**level*q0

print("CKKS in Python")
print()
print("Parameters:")
print("N: ", int(M/2), ",    q0: ", q0, ",   Scale: ", scale, ",   Levels: ", level )
z0 = [1., 1.5]
z1 = [2., 2.]

def realRounding(c):
    res = [0]*(len(c))
    for i in range(len(c)):
        res[i] = np.round(float(mp.re(c[i])), roundCS)
    return res

st = time.time()
encoder = scheme.CKKSScheme(M=M, h=h, scale=scale, P=P, level=level, q0=q0, arithmetic64=arithmetic64)
secret = encoder.keyGen(sparce_ternary)
print("Vectors:")
print("Vector 1: ", z0)
print("Vector 2: ", z1)
print()
print("Approx Bits:", np.log2(z1[1]*scale))
print()
m0, nu0 = encoder.encode(z0, scale)
m1, nu1 = encoder.encode(z1, scale)
c0 = encoder.encrypt(m0, nu0)
c1 = encoder.encrypt(m1, nu1)

print("Resultado real vs resultado CKKS")
print()
print("Suma: c1+c2")
cadd = encoder.Cadd(c0,c1)
madd = encoder.decrypt(cadd, secret)
zadd = encoder.decode(madd, scale)

zaddreal = realRounding(zadd)
print(z0[0]+z1[0], "=?", zaddreal[0], "   ,   ", z0[1]+z1[1], "=?", zaddreal[1])

print()
print("Multiplicacion por textoplano: c1*m2")
cmul_cte = encoder.Cmul_cte(c0, m1)
cmul_cte = encoder.reScale(cmul_cte)
mmul_cte = encoder.decrypt(cmul_cte, secret)
zmul_cte = encoder.decode(mmul_cte, scale)
zmulReal_cte = realRounding(zmul_cte)
z0z1 = [mp.nint(z0[0]*z1[0]), mp.nint(z0[1]*z1[1])]
print(z0z1[0], "=?", zmulReal_cte[0], "   ,   ", z0z1[1], "=?", zmulReal_cte[1])

print()
print("Multiplicacion: c1*c2")
cmul = encoder.Cmul(c0, c1)
cmulRescale = encoder.reScale(cmul)
mmulRescale = encoder.decrypt(cmulRescale, secret)
zmulRescale = encoder.decode(mmulRescale, scale)
zmulReScaleReal = realRounding(zmulRescale)
print(z0z1[0], "=?", zmulReScaleReal[0], "   ,   ", z0z1[1], "=?", zmulReScaleReal[1])
print()

print("Suma con nivelacion de level: c1+(c1*c2)")
cadd2 = encoder.Cadd(c0, cmulRescale)
maddRescale2 = encoder.decrypt(cadd2, secret)
zaddRescale2 = encoder.decode(maddRescale2, scale)
zaddReScaleReal2 = realRounding(zaddRescale2)
print(z0[0]+z0z1[0], "=?", zaddReScaleReal2[0], "   ,   ", z0[1]+z0z1[1], "=?", zaddReScaleReal2[1])


print()
print("Multiplicacion con nivelacion de level: c1*(c1*c2)")
cmul2 = encoder.Cmul(c0, cmulRescale)
cmulRescale2 = encoder.reScale(cmul2)
#print("ql", scale*q0)
#print(cmulRescale2)
mmulRescale2 = encoder.decrypt(cmulRescale2, secret)
#print(mmulRescale2)
zmulRescale2 = encoder.decode(mmulRescale2, scale)
zmulReScaleReal2 = realRounding(zmulRescale2)
print(z0[0]*z0z1[0], "=?", zmulReScaleReal2[0], "   ,   ", z0[1]*z0z1[1], "=?", zmulReScaleReal2[1])

print()
print("Multiplicacion Otro nivel: (c1*c2)*(c1*c2)")
cmul3 = encoder.Cmul(cmulRescale, cmulRescale)
cmulRescale3 = encoder.reScale(cmul3)
print(cmulRescale3.B)
mmulRescale3 = encoder.decrypt(cmulRescale3, secret)
zmulRescale3 = encoder.decode(mmulRescale3, scale)
zmulReScaleReal3 = realRounding(zmulRescale3)
print(z0z1[0]*z0z1[0], "=?", zmulReScaleReal3[0], "   ,   ", z0z1[1]*z0z1[1], "=?", zmulReScaleReal3[1])


et = time.time()
elapsed_time = round(et - st,5)
print('Execution time:', elapsed_time, 'seconds')
