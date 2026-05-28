import lab2
import numpy as np
import matplotlib.pyplot as plt


def source(x: float) -> float:
    return (2 * np.pi)**2 * np.sin(2 * np.pi * x)


# def u_true(x):
#     return - 2 * np.pi * np.cos(2 * np.pi * x)


# def p_true(x):
#     return np.sin(2 * np.pi * x)


# def compute_rate(error, N_list):
#     return np.log(error[:-1] / error[1:]) / np.log(N_list[1:] / N_list[:-1])


N = 100
stokes = lab2.Stokes()

grid = lab2.Grid(N)
u_sol, p_sol = stokes.solve_problem(grid, source)

pwconstant_x = np.repeat(grid.x, 2)[1:-1]
pwconstant_p = np.repeat(p_sol, 2)

plt.plot(grid.x, u_sol)
plt.plot(pwconstant_x, pwconstant_p)
plt.legend(["Flux", "Pressure"])
# plt.show()

# plt.spy(spp)
plt.savefig('plot.png', dpi=150, bbox_inches='tight')

plt.close()

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
