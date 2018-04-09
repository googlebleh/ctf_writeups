#!/usr/bin/python

import binascii
import itertools

import gmpy2


# parse challenge
e = 3  # from problem description

with open("moduli.txt") as f:
    moduli_strs = eval(f.read())
    moduli = map(lambda s: int(s, base=16), moduli_strs)

with open("ciphertexts.txt") as f:
    ct_strs = eval(f.read())
    cts = map(lambda s: int(s, base=16), ct_strs)


# CRT implementation from
#   https://rosettacode.org/wiki/Chinese_remainder_theorem#Python
# made more efficient with gmpy2
def chinese_remainder(n, a):
    sum = 0
    prod = reduce(gmpy2.mul, n)
 
    for n_i, a_i in zip(n, a):
        p = gmpy2.div(prod, n_i)
        sum += gmpy2.mul(gmpy2.mul(a_i, mul_inv(p, n_i)), p)
    return sum % prod

def mul_inv(a, b):
    b0 = b
    x0, x1 = 0, 1
    if b == 1: return 1
    while a > 1:
        q = gmpy2.div(a, b)
        a, b = b, a%b
        x0, x1 = x1 - gmpy2.mul(q, x0), x0
    if x1 < 0:
        x1 = gmpy2.add(x1, b0)
    return x1


def main():
    # make sure moduli are pairwise relatively prime
    for i in xrange(len(moduli)):
        for j in xrange(i+1, len(moduli)):
            assert(gmpy2.gcd(moduli[i], moduli[j]) == 1)

    # decrypt using Hastad's attack
    # https://cims.nyu.edu/~regev/teaching/lattices_fall_2004/ln/rsa.pdf
    for ct_triple in itertools.combinations(cts, 3):
        pt_cubed = chinese_remainder(moduli, ct_triple)
        pt, cuberoot_ok = gmpy2.iroot(pt_cubed, 3)
        if cuberoot_ok:
            print binascii.unhexlify("{:x}".format(pt))
            return 0

    return 1

 
if __name__ == '__main__':
    main()
