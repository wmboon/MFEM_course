import ngsolve as ng
from netgen.occ import unit_cube
import numpy as np


def test_preconditioner(K_list, h, omega=1):

    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=h))

    V = ng.HDiv(mesh, order=0, RT=True)
    Q = ng.L2(mesh, order=0)
    fes = V * Q

    (u, p), (v, q) = fes.TnT()

    K = 1.0

    b = ng.LinearForm(fes)
    # normal = ng.specialcf.normal(3)
    b += v.Trace() * ng.x * ng.ds
    b.Assemble()

    iter_list = []

    for K in K_list:
        print(f"    {K=:.2e}")
        K_inv = 1 / K

        a = ng.BilinearForm(fes)
        a += K_inv * u * v * ng.dx
        a += -p * ng.div(v) * ng.dx
        a += -q * ng.div(u) * ng.dx

        a.Assemble()

        precond = ng.BilinearForm(fes)
        precond += K_inv * u * v * ng.dx
        precond += (omega + K_inv) * ng.div(u) * ng.div(v) * ng.dx
        precond += 1 / (omega + K_inv) * p * q * ng.dx

        # p_multigrid = ng.Preconditioner(precond, "multigrid")
        precond.Assemble()

        P = precond.mat.Inverse()

        solver = ng.solvers.MinResSolver(mat=a.mat, pre=P)

        _ = solver.Solve(b.vec)

        iter_list.append(solver.iterations)

    return iter_list


if __name__ == "__main__":
    h_list = 0.5 ** np.arange(1, 4)
    K_list = 10.0 ** np.arange(-6, 7, 2)

    iter_table = []

    for h in h_list:
        print(f"{h=:.2e}")
        iter_table.append(test_preconditioner(K_list, h, 0))

    iter_array = np.array(iter_table)
    print(iter_array)

    pass
