import numpy as np
import matplotlib.pyplot as plt

# Copied from lab1


def solve_FEM_Poisson_non_homogeneous(N):
    h = 1/N
    x = np.linspace(0, 1, N + 1)
    A_E = np.array([[1, -1], [-1, 1]]) / h

    A = np.zeros((N + 1, N + 1))

    for i in range(N):
        loc_dofs = [i, i+1]
        A[np.ix_(loc_dofs, loc_dofs)] += A_E

    b = np.zeros(N + 1)

    for i in range(N):
        loc_dofs = [i, i+1]
        b[loc_dofs] += h / 2

    u = np.zeros(N + 1)
    u[0] = 1

    b -= A @ u

    freedofs = np.arange(1, N + 1)
    A_free = A[np.ix_(freedofs, freedofs)]
    b_free = b[freedofs]

    b_free[-1] += -1

    u_free = np.linalg.solve(A_free, b_free)

    u[freedofs] = u_free

    return u, x


u, x = solve_FEM_Poisson_non_homogeneous(10)

plt.plot(x, u)
plt.show()


# Mass matrix computation
N = 10
h = 1/N

M_E = np.array([[2, 1], [1, 2]]) * h / 6

M = np.zeros((N + 1, N + 1))

for i in range(N):
    loc_dofs = [i, i+1]
    M[np.ix_(loc_dofs, loc_dofs)] += M_E

print(M * 6)
