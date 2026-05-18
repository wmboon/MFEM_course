import lab2
import numpy as np
import matplotlib.pyplot as plt


def source(x):
    return (2 * np.pi) ** 2 * np.sin(2 * np.pi * x)


def compute_rate(error, N_list):
    return np.log(error[:-1] / error[1:]) / np.log(N_list[1:] / N_list[:-1])


mixeddarcy = lab2.MixedDarcyRobin()


N = 300
grid = lab2.Grid(N)
u_sol, p_sol = mixeddarcy.solve_problem(grid, source)


pwconstant_x = np.repeat(grid.x, 2)[1:-1]
pwconstant_p = np.repeat(p_sol, 2)

plt.plot(grid.x, u_sol)
plt.plot(pwconstant_x, pwconstant_p)

plt.show()
