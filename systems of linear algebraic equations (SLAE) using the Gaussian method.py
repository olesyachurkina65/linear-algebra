import numpy as np
from fractions import Fraction
from numpy.linalg import norm

def gauss_float(A, b):
    A, b = A.copy(), b.copy()
    n = len(b)
    for i in range(n):
        for j in range(i+1, n):
            f = A[j, i] / A[i, i]
            A[j, i:] -= f * A[i, i:]
            b[j] -= f * b[i]
    x = np.zeros(n)
    for i in reversed(range(n)):
        x[i] = (b[i] - A[i, i+1:] @ x[i+1:]) / A[i, i]
    return x

def gauss_fraction(A, b):
    A = [row.copy() for row in A]
    b = b.copy()
    n = len(b)
    for i in range(n):
        for j in range(i+1, n):
            f = A[j][i] / A[i][i]
            A[j][i:] = [A[j][k] - f * A[i][k] for k in range(i, n)]
            b[j] -= f * b[i]
    x = [Fraction(0) for _ in range(n)]
    for i in reversed(range(n)):
        s = sum(A[i][j] * x[j] for j in range(i+1, n))
        x[i] = (b[i] - s) / A[i][i]
    return x

def hilbert(n, use_fraction=False):
    if use_fraction:
        return [[Fraction(1, i + j + 1) for j in range(n)] for i in range(n)]
    return np.array([[1 / (i + j + 1) for j in range(n)] for i in range(n)])

def experiment(n):
    H = hilbert(n)
    b = np.ones(n)
    x_f = gauss_float(H, b)

    Hf = hilbert(n, use_fraction=True)
    bf = [Fraction(1) for _ in range(n)]
    x_e = gauss_fraction(Hf, bf)
    x_e_f = np.array([float(xi) for xi in x_e])

    r = H @ x_e_f - b
    e = x_f - x_e_f

    print(f"n = {n}:")
    print(f"   Норма нев'язки: {norm(r):.2e}")
    print(f"   Норма різниці:  {norm(e):.2e}")
    print("")

if __name__ == "__main__":
    for n in range(10, 13):
        experiment(n)