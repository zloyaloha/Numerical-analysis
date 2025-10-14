from math import sqrt


def dot_to_matrix(v):
    res = []
    for i in range(len(v)):
        string = []
        for j in range(len(v)):
            string += [v[i] * v[j]]
        res += [string]
    return res


def dot_to_scalar(v):
    ans = 0
    for i in range(len(v)):
        ans += (v[i] * v[i])
    return ans


def delta(m1, m2):
    res = []
    for i in range(len(m1)):
        string = []
        for j in range(len(m1[i])):
            string += [m1[i][j] - m2[i][j]]
        res += [string]
    return res


def dot_on_num(num, matrix):
    res = []
    for i in range(len(matrix)):
        string = []
        for j in range(len(matrix[i])):
            string += [num * matrix[i][j]]
        res += [string]
    return res


def dot_matrix_matrix(n, m1, m2):
    res = []
    for i in range(n):
        string = []
        for j in range(n):
            ans = 0
            for k in range(n):
                ans += (m1[i][k] * m2[k][j])
            string += [ans]
        res += [string]
    return res


n = int(input())  # размер квадратной матрицы
a = [[float(j) for j in input().split()] for i in range(n)]  # квадратная матрица
E = float(input())  # точность, с которой должны быть найдены собственные значения

a_now = []  # матрица a, с которой будет происходить работа на каждом шаге, в итоге будет диагональная матрица
for i in a:
    a_now += [i[:]]

delta1 = 0  # сходимость алгоритма
delta2 = 0  # наличие мнимых собственных значений
for i in range(n):
    for j in range(i + 1, n):
        if j != i + 1:
            delta1 += (a_now[j][i] ** 2)
        delta2 += (a_now[j][i] ** 2)

while sqrt(delta1) >= E:
    # строим QR-разложение
    Q = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    R = []
    for i in range(n - 1):
        H = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
        v = []
        for j in range(i):
            v += [0]
        sign_num = a_now[i][i] / abs(a_now[i][i])
        norma = 0
        for j in range(i, n):
            norma += (a_now[j][i] ** 2)
        v += [a_now[i][i] + sign_num * (norma ** 0.5)]
        for j in range(i + 1, n):
            v += [a_now[j][i]]
        # print(v)
        H = delta(H, dot_on_num(2 / dot_to_scalar(v), dot_to_matrix(v)))
        a_now = dot_matrix_matrix(n, H, a_now)
        Q = dot_matrix_matrix(n, Q, H)
        print(H)
        print(a_now)
        print(Q)
        print("---")

    for i in a_now:
        R += [i[:]]

    a_now = dot_matrix_matrix(n, R, Q)

    # print(a_now)

    delta1 = 0
    delta2 = 0
    for i in range(n):
        for j in range(i + 1, n):
            if j != i + 1:
                delta1 += (a_now[j][i] ** 2)
            delta2 += (a_now[j][i] ** 2)

    # print(delta1)
    # break
