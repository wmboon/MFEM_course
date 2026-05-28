import lab2
import numpy as np
import matplotlib.pyplot as plt


def source(x):
    return (2 * np.pi) ** 2 * np.sin(2 * np.pi * x)


def compute_rate(error, N_list):
    return np.log(error[:-1] / error[1:]) / np.log(N_list[1:] / N_list[:-1])


N = 300
grid = lab2.Grid(N)
stokes = lab2.Stokes()
u_stokes, p_stokes = stokes.solve_problem(grid, source)


labda_list = 10 ** np.arange(0, 10)

H1_norm = stokes.assemble_mass_matrix(grid)
error_list = np.array([])

for labda in labda_list:
    herrmann = lab2.Herrmann_elasticity(labda=labda)
    u_sol, p_sol = herrmann.solve_problem(grid, source)
    diff = u_sol - u_stokes
    error_u = np.sqrt(diff.T @ H1_norm @ diff)
    error_list = np.append(error_list, error_u)

print(error_list)
pwconstant_x = np.repeat(grid.x, 2)[1:-1]
pwconstant_p = np.repeat(p_sol, 2)

plt.loglog(labda_list, error_list)

# plt.plot(grid.x, u_sol)
# plt.plot(pwconstant_x, pwconstant_p)

plt.show()
