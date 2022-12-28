import numpy as np
import mpmath as mp
# To print nicely the starting modulus qL
import locale
locale.setlocale(locale.LC_ALL, '')
from dataclasses import dataclass


"""
ql = p^l * q0, solo para  0<l!!!
Paper use: N >= (lambda+110)/7.2 * log(P*qL)
P is similar to qL.
Sigma 3.2 y h=64

Parameters that works qo=2**4, scale = 2**9 and P = scale**level*q0.
With M=8, h=2 and sigma = 3.2
"""

sigma = 3.2
LEVELING_RESCALE = False

@dataclass
class cyphertext:
    first: list
    second: list
    level: int
    nu: int
    B: int


class CKKSScheme:
    """
    Basic CKKS implementation without Bootstrapping and RNS.
    Also the rotation of ciphertexts is not implemented.
    """

    def __init__(self, M:int, h:int, scale:int, P:int, level:int, q0:int, arithmetic64=False):
        if arithmetic64:
           import arithmetic64 as arithmetic
           import distributions64 as dist
        else:
           import arithmetic as arithmetic
           import distributions as dist

        self.arithmetic = arithmetic
        self.dist = dist
        self.xi = np.exp(2 * np.pi * 1j / M)
        self.M = M
        self.N = int(M/2)
        self.original_len = int(self.N/2)
        self.h = h
        self.poly_mod = np.array([1] + [0] * (self.N - 1) + [1])
        self.level_max = level
        self.q0 = q0
        self.scale = scale
        self.sigma = sigma
        self.create_sigma_R_basis()
        self.P = P
        self.pub = ()
        self.evk = ()
        print("Starting modulus qL: ", "{0:n}".format(self.new_mod(self.level_max)))
        print()
        qL = self.new_mod(level)
        cota_N = int(110/7.2*mp.log(P*qL, 2))
        print("Security condition:", self.N, ">=", cota_N)
        self.B_clean = 8 *np.sqrt(2)*sigma*self.N + 6*sigma*np.sqrt(self.N) + 16*sigma*np.sqrt(h*self.N)
        print("scale: ", scale, " >? ",int( self.N + 2*self.B_clean))
        self.B_sk = 8*sigma*self.N/np.sqrt(3)
        self.B_scale = np.sqrt(self.N/3)*(3+8*np.sqrt(h))

    def new_mod(self, level):
        ql = (self.scale**level)*self.q0
        return ql


    @staticmethod
    def vandermonde_matrix(xi: np.complex128, M: int) -> list:
        """Computes the Vandermonde matrix from a m-th root of unity."""
        N = M //2
        matrix = []
        # We will generate each row of the matrix
        for i in range(N):
            # For each row we select a different root
            root = xi ** (2 * i + 1)
            row = []
            # Then we store its powers
            for j in range(N):
                row.append(root ** j)
            matrix.append(row)
        return matrix


    def sigma_inverse(self, b: list):
        """Encodes the vector b in a polynomial using an M-th root of unity."""
        A = self.vandermonde_matrix(self.xi, self.M)
        coef = np.linalg.solve(A, b)
        coef = self.arithmetic.rint(np.round(np.real(coef)))
        return coef


    def _sigma(self, coef) -> list:
        """Decodes a polynomial by applying it to the M-th roots of unity."""
        outputs = self.arithmetic.roots(coef, self.N, self.xi)

        return outputs


    def pi_projection(self, z: list) -> list:
        """Projects a vector of H into C^{N/2}."""
        res = z[:self.N]
        return res


    def pi_projection_inverse(self, z: list) -> list:
        """Expands a vector of C^{N/2} by expanding it with its
        complex conjugate."""

        z_conjugate = z[::-1]
        z_conjugate = [np.conjugate(x) for x in z_conjugate]
        return np.concatenate([z, z_conjugate])


    def create_sigma_R_basis(self):
        """Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))."""

        self.sigma_R_basis = np.array(self.vandermonde_matrix(self.xi, self.M)).T


    def compute_basis_coordinates(self, z):
        """Computes the coordinates of a vector with respect to the orthogonal lattice basis."""
        output = np.array([np.real(np.vdot(z, b) / np.vdot(b,b)) for b in self.sigma_R_basis])
        return output


    def sigma_R_discretization(self, z):
        """Projects a vector on the lattice using coordinate wise random rounding."""
        coordinates = self.compute_basis_coordinates(z)

        rounded_coordinates = self.arithmetic.coordinate_wise_random_rounding(coordinates)
        y = np.matmul(self.sigma_R_basis.T, rounded_coordinates)
        return y


    def pad_zeros(self, z:list):
        """ Expand vector to right size.
        It only saves the original size of one vector...
        """
        self.original_len = len(z)
        while(len(z)!=int(self.N/2)):
            z.append(0)

        return z

    def encode(self, z: list, scale):
        """Encodes a vector by expanding it first to H,
        scale it, project it on the lattice of sigma(R), and performs
        sigma inverse.
        """
        if(len(z)!=int(self.N/2)):
            z = self.pad_zeros(z)
        pi_z = self.pi_projection_inverse(z)
        scaled_pi_z = self.arithmetic.mul_cte(pi_z, scale)
        rounded_scale_pi_zi = self.sigma_R_discretization(scaled_pi_z)
        coef = self.sigma_inverse(rounded_scale_pi_zi)
        coef = self.arithmetic.rint(np.real(coef))
        nu =np.linalg.norm(coef)
        return (coef, nu)


    def decode(self, coef, scale) -> list:
        """decodes a polynomial by removing the scale,
        evaluating on the roots, and project it on c^(n/2)"""
        rescaled_p = self.arithmetic.div_list(coef, scale)
        z = self._sigma(rescaled_p)
        pi_z = self.pi_projection(z)
        pi_z_unpadding = pi_z[:self.original_len]
        return pi_z_unpadding


    def _pub_keygen(self, ql, s, poly_mod):
        """ generates the public key."""
        a = self.dist.sampleRing(ql, self.N)
        e = self.dist.sampleDG(self.sigma, self.N)
        neg_a = self.arithmetic.mul_cte(a, -1)
        axs = self.arithmetic.polyMul(neg_a, s, ql, poly_mod)
        b = self.arithmetic.polyAdd(axs, e, ql)
        pub = (b, a)
        return pub


    def _evk_keygen(self,ql ,s, poly_mod):
        """ generates the evaluation key."""
        p = self.P
        pql = p*ql
        aprim = self.dist.sampleRing(pql, self.N)
        eprim = self.dist.sampleDG(self.sigma, self.N)
        neg_aprim = self.arithmetic.mul_cte(aprim, -1)
        axsprim = self.arithmetic.polyMul(neg_aprim, s, pql, poly_mod)
        ps = self.arithmetic.mul_cte(s, p)
        ps2 = self.arithmetic.polyMul(ps, s, pql, poly_mod)

        bprim = self.arithmetic.polyAdd(axsprim, eprim, pql)
        bprim = self.arithmetic.polyAdd(bprim, ps2, pql)
        evk = (bprim, aprim)
        return evk


    def keyGen(self, sparce_ternary=False):#, _lambda=0):
        """ generates the public, secret and evaluation key.
        the security parameter lambda is not implemented.
        """
        ql = self.new_mod(self.level_max)
        poly_mod = self.poly_mod

        # secret key
        s = self.dist.sampleSecret(self.h, self.N, sparce_ternary=sparce_ternary)
        secret = (1, s)

        # pub key
        self.pub = self._pub_keygen(ql=ql, s=s, poly_mod=poly_mod)

        # evaluation key
        self.evk = self._evk_keygen(ql=ql, s=s, poly_mod=poly_mod)
        return secret


    def encrypt(self, m, nu):
        ql = self.new_mod(self.level_max)
        poly_mod = self.poly_mod
        pub = self.pub
        e0 = self.dist.sampleDG(self.sigma, self.N)
        e1 = self.dist.sampleDG(self.sigma, self.N)

        v = self.dist.sampleZO(self.N, 0.5)
        mc = (self.arithmetic.polyAdd(m, e0, ql), e1)

        vpk0 = self.arithmetic.polyMul(v, pub[0], ql, poly_mod)
        vpk1 = self.arithmetic.polyMul(v, pub[1], ql, poly_mod)

        c0 = self.arithmetic.polyAdd(vpk0, mc[0], ql)
        c1 = self.arithmetic.polyAdd(vpk1, mc[1], ql)
        c = cyphertext(c0, c1, self.level_max, nu, self.B_clean)
        return c


    def decrypt(self, c: cyphertext, secret):
        level = c.level
        ql = self.new_mod(level)
        poly_mod = self.poly_mod

        s = secret[1]
        b = c.first
        a = c.second

        axs = self.arithmetic.polyMul(a, s, ql, poly_mod, center_mod=True)
        mu = self.arithmetic.polyAdd(b, axs, ql, center_mod=True)
        return mu

    def cypher_leveling_mod(self, c:cyphertext, level):
        ql2 = self.new_mod(level)
        c0 = self.arithmetic.mod(c.first, ql2)
        c1 = self.arithmetic.mod(c.second, ql2)
        cl = cyphertext(c0, c1, level, c.nu, c.B)
        return cl


    def leveling_rescale(self, c1:cyphertext, c2:cyphertext):
        """ rescale the ciphertexts with great level depth."""
        if c1.level>c2.level:
            while(c1.level!=c2.level):
                c1 = self.reScale(c1)
        else:
            while(c1.level!=c2.level):
                c2 = self.reScale(c2)
        return c1, c2


    def leveling_mod(self, c1:cyphertext, c2:cyphertext):
        """ take modulus of the ciphertexts with lower level depth,  to
        the ciphertexts with greater level depth.
        """
        if c1.level>c2.level:
            while(c1.level!=c2.level):
                c1 = self.cypher_leveling_mod(c1, c2.level)
        else:
            while(c1.level!=c2.level):
                c2 = self.cypher_leveling_mod(c2, c1.level)
        return c1, c2


    def leveling(self, c1:cyphertext, c2:cyphertext, leveling_rescale=False):
        """ the paper gives two versions of how to operate with two ciphertexts
        that haves different level depth.
        in the paper they use simple modulus reduction.
        """
        if leveling_rescale:
            c1, c2 = self.leveling_rescale(c1, c2)
        else:
            c1, c2 = self.leveling_mod(c1, c2)
        return c1, c2


    def Cadd(self, c1: cyphertext, c2: cyphertext):
        """ add two ciphertexts of the same level.
        if they are in different level the one with greater level is level down.
        """
        if c1.level!=c2.level:
            c1, c2 = self.leveling(c1, c2, LEVELING_RESCALE)

        level = c1.level
        ql = self.new_mod(level)
        b = self.arithmetic.polyAdd(c1.first, c2.first, ql)
        a = self.arithmetic.polyAdd(c1.second, c2.second, ql)
        c = cyphertext(b, a, level, c1.nu + c2.nu, c1.B+c2.B)
        return c


    def Cmul_cte(self, c: cyphertext, m: tuple):
        """ multiply one ciphertexts two a plaintext."""
        level = c.level
        ql = self.new_mod(level)
        poly_mod = self.poly_mod
        b = self.arithmetic.polyMul(m, c.first, ql, poly_mod)
        a = self.arithmetic.polyMul(m, c.second, ql, poly_mod)
        cmult = cyphertext(b, a, level, c.nu, c.B)
        return cmult


    def relinearization(self, cmult: tuple):
        """ after a multiplication the size of the ciphertexts grows, and for
        decrypted its need the square of the secret key.
        the relinealization transform the ciphertexts of a multiplication to
        the original size.
        """
        evk = self.evk
        level = cmult[3]
        ql = self.new_mod(level)
        poly_mod = self.poly_mod

        d0 = cmult[0]
        d1 = cmult[1]
        d2 = cmult[2]
        d2p = self.arithmetic.div_list(d2,self.P)
        level = cmult[3]
        p0 = self.arithmetic.polyMul(d2p, evk[0], ql, poly_mod)
        p1 = self.arithmetic.polyMul(d2p, evk[1], ql, poly_mod)
        rel0 = self.arithmetic.polyAdd(d0, p0, ql)
        rel1 = self.arithmetic.polyAdd(d1, p1, ql)

        cmultr = (rel0, rel1, level)
        return cmultr


    def Cmul(self, c1: cyphertext, c2: cyphertext):
        """ multiply two ciphertexts of the same level.
        if they are in different level the one with greater level is level down.
        """
        if c1.level!=c2.level:
            c1, c2 = self.leveling(c1, c2, LEVELING_RESCALE)

        level = c1.level
        ql = self.new_mod(level)
        poly_mod = self.poly_mod
        b1 = c1.first
        a1 = c1.second
        b2 = c2.first
        a2 = c2.second
        d0 = self.arithmetic.polyMul(b1, b2, ql, poly_mod)
        a1b2 = self.arithmetic.polyMul(a1, b2, ql, poly_mod)
        a2b1 = self.arithmetic.polyMul(a2, b1, ql, poly_mod)
        d1 = self.arithmetic.polyAdd(a1b2, a2b1, ql)
        d2 = self.arithmetic.polyMul(a1, a2, ql, poly_mod)
        cmult = (d0, d1, d2, level)
        cmultr = self.relinearization(cmult)

        B_mul = (ql*self.B_sk)/self.P+self.B_scale
        B_total = c1.nu*c2.B + c1.B*c2.nu + c1.B*c2.B + B_mul
        cmultR = cyphertext(cmultr[0], cmultr[1], level, c1.nu*c2.nu, B_total)
        return cmultR


    def reScale(self, c: cyphertext):
        """ when two ciphertexts are multiply, the scaling factor of the two
        also multiplys, making a big grow in the ciphertext and error.
        the rescaling transform the result keeping a linear scaling factor.
        """
        level = c.level
        ql1 = self.new_mod(level)
        level -=  1
        ql_1 = self.new_mod(level)
        ql_1_div_ql = self.arithmetic.div(ql_1,ql1)
        c0 = self.arithmetic.mul_cte(c.first, ql_1_div_ql)
        c0 = self.arithmetic.rint(c0)
        c1 = self.arithmetic.mul_cte(c.second, ql_1_div_ql)
        c1 = self.arithmetic.rint(c1)

        B_scale = np.sqrt(self.N/3) * (3 + 8*np.sqrt(self.h))
        B = ql_1_div_ql*c.B + B_scale
        nu = ql_1_div_ql*c.nu
        c_rescaleado = cyphertext(c0, c1, level, nu, B)
        return c_rescaleado
