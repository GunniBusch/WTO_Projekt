import numpy as np
from numpy.polynomial import polynomial as P

from BFV import CT


class BFV:
    def __init__(self, ring_size: int, koef_pt: int, koef_ct: int, poly_modulo: np.ndarray,
                 modulo_relin: np.ndarray, mean, scale, encode_base):
        """

        :param ring_size: Größe des Polynomrings
        :param koef_pt: Koeffizienten des Ausgangstextes
        :param koef_ct: Koeffizienten des Chiffretextes
        :param poly_modulo: Modulo Polynom
        :param mean: Median
        :param scale: Standardabweichung
        """
        # Ring größe
        self.n = ring_size
        # Chiffretext polynom Koeffizienten
        self.q = koef_ct
        # Ausgangstext polynom koeffizienten
        self.t = koef_pt

        # (x^n+1)
        self.poly_modulo = poly_modulo
        self.p = modulo_relin
        self.mean = mean
        self.scale = scale
        self.encode_base = encode_base

    # Polynom Rechnungen
    @staticmethod
    def _rnd(x):

        return np.int64(np.round(x))

    def add(self, poly1: np.ndarray, poly2: np.ndarray) -> np.ndarray:

        xy = P.polyadd(poly1, poly2)

        return P.polydiv(xy, self.poly_modulo)[1]

    def mul(self, poly1: np.ndarray, poly2=None) -> np.ndarray:

        if poly2 is None:
            poly2 = poly1

        xy = P.polymul(poly1, poly2)

        return P.polydiv(xy, self.poly_modulo)[1]

    def _genPoly(self):
        """

        :return: Koeffizienten eines Polynoms des Grades n-1
        """
        rng = np.random.default_rng()
        return rng.integers(low=-1, high=2, size=self.n)

    def _genNormalPoly(self):
        """

        :return: Koeffizienten in einer Normalverteilung eines Polynoms des Grades n-1
        """
        rng = np.random.default_rng()
        return np.int64(rng.normal(self.mean, self.scale, self.n))

    def _genUniPoly(self, mod=None):
        """

        :return: Koeffizienten in Z_n eines Polynoms des Grades n-1
        """

        if mod is None:
            mod = self.q

        rng = np.random.default_rng(np.random.PCG64DXSM())
        return rng.integers(low=0, high=mod, size=self.n)

    def generateSecretKey(self):
        """
        Erzeugt privaten Schlüssel
        :return: privaten Schlüssel und öffentliche Schlüssel (sk,(pk1,pk2))
        """
        sk = self._genPoly()

        return sk

    def generatePublicKey(self, sk):
        """
        Erzeugt öffentlichen Schlüssel
        :param sk: Geheimer Schlüssel
        :return: privaten Schlüssel und öffentliche Schlüssel (sk,(pk1,pk2))
        """

        a = self._genUniPoly()
        e = self._genNormalPoly()
        # pk1= -1(pk2*sk+e)=(pk2*(-sk)-e)
        pk1 = self._rnd(self.add(self.mul(-a, sk), -e) % self.q)
        pk2 = a

        return pk1, pk2

    def generateRelinKey(self, sk):

        relin_key2 = self._genUniPoly(self.p * self.q)
        e = self._genNormalPoly()
        sk2 = P.polymul(sk, sk)
        secret_part = self.p * sk2

        relin_key1 = (self.add(self.add(
            self.mul(-relin_key2, sk),
            -e), secret_part)) % (self.p * self.q)

        return relin_key1, relin_key2

    def encode(self, n, b):
        """
        Encodiert eine Zahl.
        :param n: zu encodierende Zahl
        :param b: Basis
        :return: Koeffizienten der Basis b Repräsentation von n
        """

        return np.array([n] + [0] * (self.n - 1), dtype=np.int64) % self.t

    @staticmethod
    def decode(n, b):
        """
        Encodiert eine Zahl.
        :param n: zu encodierende Zahl
        :param b: Basis
        :return: Koeffizienten der Basis b Repräsentation von n
        """
        return P.polyval(b, n)

    def encrypt(self, pt, pk: tuple[np.ndarray, np.ndarray], ek: tuple[np.ndarray, np.ndarray]) -> CT:
        """

        :param ek: evaluationKey. Gleich wie Relinearisierungsschlüssel
        :param pt: Codierter Ausgangstext
        :param pk: tuple des öffentlichen SSchlüssels
        :return: tuple der Verschlüsselten texte
        """
        pk1, pk2 = pk
        pt = self.encode(pt, self.encode_base)

        u = self._genPoly()
        e1 = self._genNormalPoly()
        e2 = self._genNormalPoly()
        quotient = (self.q // self.t)

        pk1u = (self.mul(pk1, u) % self.q)
        pk2u = (self.mul(pk2, u) % self.q)

        ct1 = self._rnd(self.add(pk1u, self.add(e1, pt * quotient)) % self.q)
        ct2 = self._rnd(self.add(pk2u, e2) % self.q)

        return CT(self, (ct1, ct2), pk, ek)

    def decrypt(self, ct, sk):
        ct1, ct2 = ct.ct
        step1 = self._rnd((self.add(
            self.mul(ct2, sk), ct1
        ) % self.q))
        # Codierter Klartext
        dec = np.round(step1 / (self.q // self.t)) % self.t

        return np.int64(dec[0])
