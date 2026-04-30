import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps

# Copied from lab1


def grid(N):
    h = 1 / N
    x = np.linspace(0, 1, N + 1)
    return x, h


def assemble_mass_matrix(N, h):
    M_E = np.array([[2, 1], [1, 2]]) * h / 6

    # M = np.zeros((N + 1, N + 1))
    rows = np.array([], dtype=int)
    cols = np.array([], dtype=int)
    vals = np.array([])

    for i in range(N):
        loc_dofs = [i, i + 1]
        np.ix_(loc_dofs, loc_dofs)
        rows = np.append(rows, np.repeat(loc_dofs, 2))
        cols = np.append(cols, np.tile(loc_dofs, 2))
        vals = np.append(vals, M_E.flatten())

    return sps.csc_array((vals, (rows, cols)))


def assemble_div_matrix(N):
    B_E = np.array([-1, 1])

    rows = np.array([], dtype=int)
    cols = np.array([], dtype=int)
    vals = np.array([])

    for i in range(N):
        loc_dofs = [i, i + 1]
        rows = np.append(rows, np.repeat(i, 2))
        cols = np.append(cols, loc_dofs)
        vals = np.append(vals, B_E)

    return sps.csc_array((vals, (rows, cols)))


def assemble_SPP(N, h):
    A = assemble_mass_matrix(N, h)
    B = assemble_div_matrix(N)

    spp = sps.block_array([[A, -B.T], [B, None]]).tocsc()
    return spp


def assemble_rhs(N, x, h):
    cell_centers = (x[:-1] + x[1:]) / 2
    rhs_p = np.array([source(x) * h for x in cell_centers])
    rhs_v = np.zeros(N + 1)
    return np.hstack((rhs_v, rhs_p))


def source(x):
    return (2 * np.pi) ** 2 * np.sin(2 * np.pi * x)


def solve_problem(N):
    x, h = grid(N)

    spp = assemble_SPP(N, h)
    rhs = assemble_rhs(N, x, h)

    sol = sps.linalg.spsolve(spp, rhs)

    u_sol = sol[: N + 1]
    p_sol = sol[N + 1 :]
    return u_sol, p_sol


def u_true(x):
    return -2 * np.pi * np.cos(2 * np.pi * x)


def p_true(x):
    return np.sin(2 * np.pi * x)


def compute_error_flux(u_sol, u_true, N):
    x, _ = grid(N)
    u_interp = np.array([u_true(x_i) for x_i in x])
    diff = u_sol - u_interp
    M = assemble_mass_matrix(N, 1 / N)
    error_squared = diff @ M @ diff
    norm_squared = u_interp @ M @ u_interp
    return np.sqrt(error_squared / norm_squared)


def compute_error_pressure(p_sol, p_true, N):
    x, h = grid(N)
    cell_centers = (x[:-1] + x[1:]) / 2

    p_interp = np.array([p_true(x_i) for x_i in cell_centers])
    diff = p_sol - p_interp

    M = h * sps.eye_array(N)

    error_squared = diff @ M @ diff
    norm_squared = p_interp @ M @ p_interp
    return np.sqrt(error_squared / norm_squared)


def compute_rate(error, N_list):
    return np.log(error[:-1] / error[1:]) / np.log(N_list[1:] / N_list[:-1])


if __name__ == "__main__":
    N_list = np.array([5, 10, 20, 40, 60, 80, 100])
    error_u = np.array([])
    error_p = np.array([])
    for N in N_list:
        u_sol, p_sol = solve_problem(N)
        error_u = np.append(error_u, compute_error_flux(u_sol, u_true, N))
        error_p = np.append(error_p, compute_error_pressure(p_sol, p_true, N))
    plt.plot(N_list, error_u)
    plt.plot(N_list, error_p)
    plt.show()

    rates_p = compute_rate(error_p, N_list)
    rates_u = compute_rate(error_u, N_list)
    print(f"rates_p =", rates_p, "rates_u =", rates_u)

    # N = 40
    # x, h = grid(N)
    # pwconstant_x = np.repeat(x, 2)[1:-1]
    # pwconstant_p = np.repeat(p_sol, 2)

    # plt.plot(x, u_sol)
    # plt.plot(pwconstant_x, pwconstant_p)

    # plt.spy(spp)
    # plt.show()

    pass

    # h = 1 / N

    # M_E = np.array([[2, 1], [1, 2]]) * h / 6

    # M = np.zeros((N + 1, N + 1))

    # for i in range(N):
    #     loc_dofs = [i, i + 1]
    #     M[np.ix_(loc_dofs, loc_dofs)] += M_E

    # print(M * 6)
