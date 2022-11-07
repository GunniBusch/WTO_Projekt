import numpy as np

import BFV


class CT:
    def __init__(self, bfv: BFV, ct: tuple[np.ndarray, np.ndarray], pk, ek):

        self.bfv = bfv
        self.ct = ct
        self.pk = pk
        self.ek = ek

    def _relin(self, ct: tuple[np.ndarray, np.ndarray, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
        ct1, ct2, ct3 = ct
        ek1, ek2 = self.ek
        ct21 = np.int64(np.round(self.bfv.mul(ek1, ct3) / self.bfv.p)) % self.bfv.q
        ct22 = np.int64(np.round(self.bfv.mul(ek2, ct3) / self.bfv.p)) % self.bfv.q

        ct1 = np.int64(self.bfv.add(ct1, ct21)) % self.bfv.q
        ct2 = np.int64(self.bfv.add(ct2, ct22)) % self.bfv.q

        return ct1, ct2

    def _rnd(self, x):
        return self.bfv._rnd(x)

    def __add__(self, other):

        ct11, ct12 = self.ct
        if isinstance(other, CT):
            ct21, ct22 = other.ct
            ct1_new = self._rnd(self.bfv.add(ct11, ct21) % self.bfv.q)
            ct2_new = self._rnd(self.bfv.add(ct12, ct22) % self.bfv.q)

            return CT(self.bfv, (ct1_new, ct2_new), self.pk, self.ek)

        elif isinstance(other, int):

            ct1, ct2 = self.ct

            size = self.bfv.n - 1

            m = np.array([other] + [0] * (size), dtype=np.int64) % self.bfv.t

            d = self.bfv.q // self.bfv.t

            scaled_m = d * m

            new_ct1 = self.bfv.add(ct1, scaled_m) % self.bfv.q

            return CT(self.bfv, (new_ct1, ct2), self.pk, self.ek)
        else:
            raise ValueError(f"can only add ciphertexts or integers, but operand was {other.__class__}")

    def __mul__(self, other):

        ct11, ct12 = self.ct
        if isinstance(other, CT):
            ct21, ct22 = other.ct
            ct1 = np.int64(
                np.round(
                    self.bfv.mul(ct11, ct21) * self.bfv.t / self.bfv.q
                )) % self.bfv.q

            ct2 = np.int64(
                np.round(
                    self.bfv.add(self.bfv.mul(ct11, ct22),
                                 self.bfv.mul(ct12, ct21)
                                 ) * self.bfv.t / self.bfv.q
                )) % self.bfv.q

            ct3 = np.int64(
                np.round(
                    self.bfv.mul(ct12, ct22) * self.bfv.t / self.bfv.q)
            ) % self.bfv.q

            return CT(self.bfv, self._relin((ct1, ct2, ct3)), self.pk, self.ek)

        elif isinstance(other, int):

            ct1, ct2 = self.ct
            size = len(self.bfv.poly_modulo)

            m = np.array([other] + [0] * (size - 2), dtype=np.int64) % self.bfv.t
            ct_new_1 = self.bfv.mul(ct1, m) % self.bfv.q
            ct_new_2 = self.bfv.mul(ct2, m) % self.bfv.q

            return CT(self.bfv, (ct_new_1, ct_new_2), self.pk, self.ek)
        else:
            raise ValueError(f"can only multiply with ciphertexts or integers, but operand was {other.__class__}")
