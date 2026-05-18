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


mixeddarcy = lab2.MixedDarcy()

N_list = np.array([5, 10, 20, 40, 60, 80, 100])
error_u = np.array([])
error_p = np.array([])
for N in N_list:
    grid = lab2.Grid(N)
    u_sol, p_sol = mixeddarcy.solve_problem(grid, source)
    error_u = np.append(error_u, mixeddarcy.compute_error_flux(u_sol, u_true, grid))
    error_p = np.append(error_p, mixeddarcy.compute_error_pressure(p_sol, p_true, grid))
plt.plot(N_list, error_u)
plt.plot(N_list, error_p)
plt.show()

rates_p = compute_rate(error_p, N_list)
rates_u = compute_rate(error_u, N_list)
# print(f"rates_p =", rates_p, "rates_u =", rates_u)


pwconstant_x = np.repeat(grid.x, 2)[1:-1]
pwconstant_p = np.repeat(p_sol, 2)

plt.plot(grid.x, u_sol)
plt.plot(pwconstant_x, pwconstant_p)

plt.show()
