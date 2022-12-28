import numpy as np
import mpmath as mp
from numpy.polynomial import polynomial as poly
print("Arbitrary precision arithmetic")
center_mod = False
mp.dps = 50
mp.pretty = False
print(mp.dps)

def round_coordinates(coordinates):
    """Gives the integral rest."""
    coordinates = coordinates - np.floor(coordinates)
    return coordinates


def coordinate_wise_random_rounding(coordinates):
    """Rounds coordinates randonmly."""
    r = round_coordinates(coordinates)
    f = np.array([np.random.choice([c, c-1], 1, p=[1-c, c]) for c in r]).reshape(-1)

    rounded_coordinates = coordinates - f
    rounded_coordinates = [int(coeff) for coeff in rounded_coordinates]
    return rounded_coordinates


def roots(coef, N, xi):
    outputs = []
    # We simply apply the polynomial on the roots
    coef.reverse()

    for i in range(N):
        root = xi ** (2 * i + 1)
        output = mp.polyval(coef, root)
        outputs.append(output)
    return outputs

def mul_cte(x,y):
    res = [0]*len(x)
    for i in range(len(x)):
        res[i] = mp.fmul(x[i],y)
    return res


def div(x,y):
    return mp.fdiv(x,y)

def div_list(x,y):
    res = [0]*len(x)
    for i in range(len(x)):
        res[i] = mp.fdiv(x[i],y)
    return res

def modCenterVal(x, ql):
    if 0 <= x and x <= mp.fsub(ql,1):
        if x > mp.fdiv(ql,2):
            return mp.fsub(x, ql)
        else:
            return x
    else:
        new_x = mp.fmod(x, ql)
        return modCenterVal(new_x, ql)


def modVal(x, ql, center_mod=center_mod):
    if center_mod:
        x = modCenterVal(x, ql)
    else:
        x = mp.fmod(x, ql)
    return x


def mod(coef: list, ql, center_mod=center_mod):
    for i in range(len(coef)):
        coef[i] = modVal(coef[i], ql, center_mod)
    return coef


def rint(coef):
    for i in range(len(coef)):
        coef[i] = mp.nint(coef[i])
    return coef


def polynomial_multiplication(P, Q):
    m = len(P)
    n = len(Q)
    result = [0]*(m+n-1)
    for i in range(m):
        for j in range(n):
            prod =  mp.fmul(P[i],Q[j])
            res_ij = result[i+j]
            result[i+j] = mp.fadd(res_ij, prod)
    return result


def polynomial_add(P, Q):
    m = len(P)
    n = len(Q)
    if m!=n:
        maxlen = max(m, n)
        minlen = min(m, n)
        for i in range(minlen, maxlen):
            if m<n:
                P = np.append(P, 0)
            elif n<m:
                Q = np.append(Q, 0)
    else:
        maxlen = m
    result = [0]*(maxlen)
    for i in range(maxlen):
        result[i]=  mp.fadd(P[i],Q[i])
    return result


def polynomial_div(coef, poly_mod):
    n = len(coef)
    mod = len(poly_mod)-1
    if (n<=mod):
        return coef
    else:
        result = [0]*(mod)
        num_fractions = int(np.ceil(n/(mod)))
        while(len(coef)<mod*num_fractions):
            coef = np.append(coef,0)
        for i in range(mod):
            sum_coef = 0
            for j in range(num_fractions):
                sum_coef = mp.fadd(sum_coef, (-1)**(j)*coef[i+(mod)*j])
            result[i] = sum_coef
        return result

def polyMul(x, y, ql, poly_mod, center_mod=center_mod):

    xXy = mod(polynomial_div(polynomial_multiplication(x, y), poly_mod), ql, center_mod)
    return xXy


def polyAdd(x, y, ql,  center_mod=center_mod):

    return mod(polynomial_add(x, y), ql, center_mod)

