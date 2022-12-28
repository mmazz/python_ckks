from numpy.polynomial import polynomial as poly
import numpy as np
print("Arithmetic 64 bits")
center_mod = False


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
    p =  poly.Polynomial(coef)
    # We simply apply the polynomial on the roots
    for i in range(N):
        root = xi ** (2 * i + 1)
        output = p(root)
        outputs.append(output)
    return outputs

def mul_cte(x, cte):
    res = [0]*(len(x))
    for i in range(len(x)):
        res[i] = x[i]*cte
    return res

def div(x,y):
    return x/y

def div_list(x,y):
    return x/y

def modCenterVal(x, ql):
    if 0 <= x and x <= ql-1:
        if x > ql/2:
            return x - ql
        else:
            return x
    else:
        new_x = x % ql
        return modCenterVal(new_x, ql)


def modVal(x, ql, center_mod=center_mod):
    if center_mod:
        x = modCenterVal(x, ql)
    else:
        x = x % ql
    return x


def mod(coef: np.array, ql, center_mod=center_mod):
    for i in range(len(coef)):
        coef[i] = modVal(coef[i], ql, center_mod)
    return coef


def rint(coef):
    coef = np.rint(coef)
    return coef


def polyMul(x, y, ql, poly_mod, center_mod=center_mod):

    xXy = mod(poly.polydiv(poly.polymul(x, y), poly_mod)[1], ql, center_mod)
    return np.int64(xXy)


def polyAdd(x, y, ql,  center_mod=center_mod):

    return np.int64(mod(poly.polyadd(x, y), ql, center_mod))


