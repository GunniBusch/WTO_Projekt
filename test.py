import numpy as np

import BFV

if __name__ == "__main__":
    n = 2 ** 3
    # ciphertext modulus
    q = 2 ** 19
    # plaintext modulus
    t = 2 ** 5
    p = q ** 2 + 1

    modulo = np.array([1] + [0] * (n - 1) + [1])

    mean, scale = 0, 3.2
    base = t

    test = BFV.BFV(n, t, q, modulo, p, mean, scale, base)
    sk = test.generateSecretKey()
    pk = test.generatePublicKey(sk)
    ek = test.generateRelinKey(sk)

    pt1 = 11
    pt2 = 2

    ct1 = test.encrypt(pt1, pk, ek)
    ct2 = test.encrypt(pt2, pk, ek)
    ct3 = test.encrypt(pt2, pk, ek)

    print(f"Entschlüsselt(ct1 + ct2) = {test.decrypt(ct1 + ct2, sk)}; pt1 + pt2 = {pt1} + {pt2}")
    print(f"Entschlüsselt(ct1 * ct2) = {test.decrypt(ct1 * ct2, sk)}; pt1 * pt2 = {pt1} * {pt2}")
    print(f"Entschlüsselt(ct1 + 9) = {test.decrypt(ct1 + 9, sk)}; pt1 + 9 = {pt1} + {9}")
    print(f"Entschlüsselt(ct1 * 8) = {test.decrypt(ct1 * 8, sk)}; pt1 * 8 = {pt1} * {8}")
