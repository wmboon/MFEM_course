import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps

# Copied from lab1


def assemble_mass_matrix(N):
    h = 1/N
    x = np.linspace(0, 1, N + 1)
    M_E = np.array([[2, 1], [1, 2]]) * h / 6

    M = np.zeros((N + 1, N + 1))

    for i in range(N):
        loc_dofs = [i, i+1]
        M[np.ix_(loc_dofs, loc_dofs)] += M_E

    return M


def assemble_div_matrix(N):
    h = 1/N
    x = np.linspace(0, 1, N + 1)
    B_E = np.array([-1, 1])

    B = np.zeros((N, N + 1))

    for i in range(N):
        loc_dofs = [i, i+1]
        B[i, loc_dofs] += B_E

    return B


def assemble_SPP(N):
    A = assemble_mass_matrix(N)
    B = assemble_div_matrix(N)

    spp = sps.block_array([[A, B.T], [B, None]])

    return spp


def new_func():
    b = np.zeros(N + 1)

    for i in range(N):
        loc_dofs = [i, i+1]
        b[loc_dofs] += h / 2

    u = np.zeros(N + 1)
    u[0] = 1

    b -= M @ u

    freedofs = np.arange(1, N + 1)
    A_free = M[np.ix_(freedofs, freedofs)]
    b_free = b[freedofs]

    b_free[-1] += -1

    u_free = np.linalg.solve(A_free, b_free)

    u[freedofs] = u_free

    return u, x


if __name__ == "__main__":
    N = 10
    spp = assemble_SPP(N)

    plt.spy(spp)
    plt.show()

    pass

    # M = assemble_mass_matrix(N)

    # plt.plot(x, u)
    # plt.show()

    # # Mass matrix computation
    # N = 10
    # h = 1/N

    # M_E = np.array([[2, 1], [1, 2]]) * h / 6

    # M = np.zeros((N + 1, N + 1))

    # for i in range(N):
    #     loc_dofs = [i, i+1]
    #     M[np.ix_(loc_dofs, loc_dofs)] += M_E

    # print(M * 6)
