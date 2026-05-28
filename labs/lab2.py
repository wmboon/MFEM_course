import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps

# Copied from lab1


class Grid:
    def __init__(self, N):
        self.N = N
        self.h = 1 / N
        self.x = np.linspace(0, 1, N + 1)
        self.cell_centers = (self.x[:-1] + self.x[1:]) / 2


class MixedFiniteElement:
    def assemble_mass_matrix(self, grid: Grid):
        M_E = np.array([[2, 1], [1, 2]]) * grid.h / 6

        rows = np.array([], dtype=int)
        cols = np.array([], dtype=int)
        vals = np.array([])

        for i in range(grid.N):
            loc_dofs = [i, i + 1]
            rows = np.append(rows, np.repeat(loc_dofs, 2))
            cols = np.append(cols, np.tile(loc_dofs, 2))
            vals = np.append(vals, M_E.flatten())

        return sps.csc_array((vals, (rows, cols)))

    def assemble_div_matrix(self, grid: Grid):
        B_E = np.array([-1, 1])

        rows = np.array([], dtype=int)
        cols = np.array([], dtype=int)
        vals = np.array([])

        for i in range(grid.N):
            loc_dofs = [i, i + 1]
            rows = np.append(rows, np.repeat(i, 2))
            cols = np.append(cols, loc_dofs)
            vals = np.append(vals, B_E)

        return sps.csc_array((vals, (rows, cols)))

    def assemble_SPP(self, grid: Grid):
        A = self.assemble_A_matrix(grid)
        B = self.assemble_div_matrix(grid)

        spp = sps.block_array(
            [
                [A, -B.T],
                [B, None],
            ]
        ).tocsc()

        return spp

    def solve_problem(self, grid, source):
        spp = self.assemble_SPP(grid)
        rhs = self.assemble_rhs(grid, source)

        sol = sps.linalg.spsolve(spp, rhs)

        u_sol = sol[: grid.N + 1]
        p_sol = sol[grid.N + 1 :]

        return u_sol, p_sol

    def compute_error_flux(self, u_sol, u_true, grid):
        u_interp = np.array([u_true(x_i) for x_i in grid.x])
        diff = u_sol - u_interp

        M = self.assemble_mass_matrix(grid)

        error_squared = diff @ M @ diff

        norm_flux_squared = u_interp @ M @ u_interp
        return np.sqrt(error_squared / norm_flux_squared)

    def compute_error_pressure(self, p_sol, p_true, grid):
        cell_centers = (grid.x[:-1] + grid.x[1:]) / 2

        p_interp = np.array([p_true(x_i) for x_i in cell_centers])
        diff = p_sol - p_interp

        M = grid.h * sps.eye_array(grid.N)

        error_squared = diff @ M @ diff

        norm_pressure_squared = p_interp @ M @ p_interp
        return np.sqrt(error_squared / norm_pressure_squared)

    def assemble_stiffness_matrix(self, grid: Grid):
        M_E = 1 / grid.h * np.array([[1, -1], [-1, 1]])

        rows = np.array([], dtype=int)
        cols = np.array([], dtype=int)
        vals = np.array([])

        for i in range(grid.N):
            loc_dofs = [i, i + 1]
            rows = np.append(rows, np.repeat(loc_dofs, 2))
            cols = np.append(cols, np.tile(loc_dofs, 2))
            vals = np.append(vals, M_E.flatten())

        return sps.csc_array((vals, (rows, cols)))

    def assemble_rhs(self, grid, source) -> np.array:
        rhs_p = np.array([source(x) * grid.h for x in grid.cell_centers])
        rhs_v = np.zeros(grid.N + 1)

        return np.hstack((rhs_v, rhs_p))


class MixedDarcy(MixedFiniteElement):
    def assemble_A_matrix(self, grid):
        return self.assemble_mass_matrix(grid)

    # def new_func():
    #     b = np.zeros(N + 1)

    #     for i in range(N):
    #         loc_dofs = [i, i + 1]
    #         b[loc_dofs] += h / 2

    #     u = np.zeros(N + 1)
    #     u[0] = 1

    #     b -= M @ u

    #     freedofs = np.arange(1, N + 1)
    #     A_free = M[np.ix_(freedofs, freedofs)]
    #     b_free = b[freedofs]

    #     b_free[-1] += -1

    #     u_free = np.linalg.solve(A_free, b_free)

    #     u[freedofs] = u_free

    #     return u


class MixeDarcyRobin(MixedDarcy):
    def assemble_mass_matrix(self, grid: Grid):
        M = super().assemble_mass_matrix(grid)

        rows = np.array([0, grid.N])
        cols = rows
        vals = np.ones(2)

        R = sps.csc_array((vals, (rows, cols)))

        return M + R


class Stokes(MixedDarcy):
    def assemble_A_matrix(self, grid):
        return self.assemble_stiffness_matrix(grid)

    def solve_problem(self, grid, source):
        spp = self.assemble_SPP(grid)
        rhs = self.assemble_rhs(grid, source)

        sol = np.zeros(spp.shape[0])
        sol[0] = 1

        rhs -= spp @ sol

        freedofs = np.arange(1, spp.shape[0])
        spp_free = spp[np.ix_(freedofs, freedofs)]
        rhs_free = rhs[freedofs]

        # sol = np.linalg.solve(spp_free, rhs_free)

        sol_reduced = sps.linalg.spsolve(spp_free, rhs_free)

        sol[freedofs] = sol_reduced

        u_sol = sol[: grid.N + 1]
        p_sol = sol[grid.N + 1 :]

        return u_sol, p_sol

    def assemble_rhs(self, grid: Grid, source) -> np.array:

        f_interp = [source(x) for x in grid.x]
        f_interp = np.array(f_interp)

        M = self.assemble_mass_matrix(grid)

        rhs_v = M @ f_interp
        rhs_p = np.zeros(grid.N)
        # np.array([source(x) * grid.h for x in grid.cell_centers])

        return np.hstack((rhs_v, rhs_p))


class Herrmann_elasticity(Stokes):
    def __init__(self, labda):
        self.labda = labda

    def assemble_SPP(self, grid):

        A = self.assemble_stiffness_matrix(grid)
        B = self.assemble_div_matrix(grid)
        C = self.assemble_C_matrix(grid, self.labda)

        spp = sps.block_array(
            [
                [A, -B.T],
                [B, C],
            ]
        ).tocsc()

        return spp

    def assemble_C_matrix(self, grid: Grid, labda):
        return 1 / labda * grid.h * sps.eye_array(grid.N)
