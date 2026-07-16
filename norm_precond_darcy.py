import ngsolve as ng
import numpy as np
from netgen.occ import unit_cube


def test_preconditioner(K_list, h):

    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=h))
    omega = 1

    V = ng.HDiv(mesh, order=0, RT=True)
    Q = ng.L2(mesh, order=0)
    fes = V * Q

    (u, p), (v, q) = fes.TnT()

    b = ng.LinearForm(fes)
    normal = ng.specialcf.normal(3)
    b += normal * v.Trace() * ng.x * ng.ds
    b.Assemble()

    iter_list = []

    for K in K_list:
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
        precond.Assemble()

        P = precond.mat.Inverse()
        P2 = ng.MultiGridPreconditioner(a)

        solver = ng.solvers.MinResSolver(mat=a.mat, pre=P2)
        solver.Solve(b.vec)
        # print(solver.iterations)

        iter_list.append(solver.iterations)
        print(f"         K={K}")
    return iter_list


if __name__ == "__main__":
    h_list = 0.5 ** np.arange(4, 5)
    K_list = 10.0 ** np.arange(-6, 7, 2)

    iter_table = []

    for h in h_list:
        iter_table.append(test_preconditioner(K_list, h))
        print(f"h={h} ist fertig")

    iter_array = np.array(iter_table)
    print(iter_array)

pass
