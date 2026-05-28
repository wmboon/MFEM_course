import lab2
import numpy as np
import matplotlib.pyplot as plt


def source(x):
    return (2 * np.pi) ** 2 * np.sin(2 * np.pi * x)


def u_true(x):
    return -2 * np.pi * np.cos(2 * np.pi * x)


def p_true(x):
    return np.sin(2 * np.pi * x)


def compute_rate(error, N_list):
    return np.log(error[:-1] / error[1:]) / np.log(N_list[1:] / N_list[:-1])


N = 100
mixedDarcy = lab2.MixeDarcyRobin()

grid = lab2.Grid(N)
u_sol, p_sol = mixedDarcy.solve_problem(grid, source)

pwconstant_x = np.repeat(grid.x, 2)[1:-1]
pwconstant_p = np.repeat(p_sol, 2)

plt.plot(grid.x, u_sol)
plt.plot(pwconstant_x, pwconstant_p)
plt.show()

# plt.spy(spp)
plt.savefig("plot.png", dpi=150, bbox_inches="tight")
plt.close()

pass
