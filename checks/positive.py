from math import lcm

def psd(matrix):
    n = len(matrix)
    assert all((len(r) == n for r in matrix))
    assert all((matrix[i][j] == matrix[j][i] for i in range(n) for j in range(i)))
    denom = 1
    for row in matrix:
        for v in row:
            denom = lcm(denom, v.denominator)
    a = [[v.numerator * (denom // v.denominator) for v in row] for row in matrix]
    previous = 1
    rank = 0
    for k in range(n):
        assert a[k][k] >= 0, 'negative diagonal/pivot'
        if a[k][k] == 0:
            assert all((a[k][j] == 0 for j in range(k, n))), 'zero diagonal with nonzero row'
            pivot = next((j for j in range(k + 1, n) if a[j][j] != 0), None)
            if pivot is None:
                assert all((a[i][j] == 0 for i in range(k, n) for j in range(k, n)))
                return rank
            a[k], a[pivot] = (a[pivot], a[k])
            for row in a:
                row[k], row[pivot] = (row[pivot], row[k])
        p = a[k][k]
        assert p > 0
        for i in range(k + 1, n):
            for j in range(i, n):
                value, rem = divmod(p * a[i][j] - a[i][k] * a[k][j], previous)
                assert rem == 0, 'fraction-free division is inexact'
                a[i][j] = a[j][i] = value
        previous = p
        rank += 1
    return rank
