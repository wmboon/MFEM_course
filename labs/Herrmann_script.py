import lab2
import numpy as np
import matplotlib.pyplot as plt


def source(x):
    return (2 * np.pi) ** 2 * np.sin(2 * np.pi * x)


N = 100
labda_list = 10 ** np.arange(0, 10)
grid = lab2.Grid(N)

stokes = lab2.Stokes()
u_stokes, p_stokes = stokes.solve_problem(grid, source)

H1_norm = stokes.assemble_mass_matrix(grid)
error_list = []

for labda in labda_list:
    herrmann = lab2.Herrmann_elasticity(labda)
    u_sol, p_sol = herrmann.solve_problem(grid, source)

    diff = u_sol - u_stokes
    error_list.append(np.sqrt(diff.T @ H1_norm @ diff))


errors = np.array(error_list)
print(errors)

plt.loglog(labda_list, errors)

pwconstant_x = np.repeat(grid.x, 2)[1:-1]
pwconstant_p = np.repeat(p_sol, 2)

# plt.plot(grid.x, u_sol)
# plt.plot(pwconstant_x, pwconstant_p)
plt.show()

# plt.spy(spp)
plt.savefig("plot.png", dpi=150, bbox_inches="tight")
plt.close()

pass
