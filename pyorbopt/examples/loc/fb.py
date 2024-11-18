import numpy as np
from pyscf import gto, scf, lo

from pyorbopt.pyorbopt import solver

mol = gto.Mole()
mol.build(
    atom="""
        O	 0.0000000	 0.0000000	 0.0000000
        H	 0.7569685	 0.0000000	-0.5858752
        H	-0.7569685	 0.0000000	-0.5858752
    """,
    basis="631g",
    symmetry="C2v",
    verbose=0,
)

# hf calculation
hf = scf.RHF(mol)
hf.conv_tol = 1.0e-10
hf.kernel()

# orbitals
orbs = [hf.mo_coeff[:, : min(mol.nelec)], hf.mo_coeff[:, max(mol.nelec) :]]

# loop over occupied and virtual subspaces
for mo_coeff in orbs:
    # number of orbitals
    norb = mo_coeff.shape[1]

    # atomic initial guess
    u = lo.boys.atomic_init_guess(mol, mo_coeff)
    mo_coeff = mo_coeff @ u

    # PySCF FB localization
    loc = lo.Boys(mol, mo_coeff)
    loc.conv_tol = 1.0e-10
    loc.verbose = 4
    loc.kernel(mo_coeff=mo_coeff)

    # initialize pyscf localization object
    loc = lo.Boys(mol, mo_coeff)

    # unpack matrix
    def unpack(kappa):
        matrix = np.zeros(2 * (norb,), dtype=np.float64)
        idx = np.tril_indices(norb, -1)
        matrix[idx] = kappa
        return matrix - matrix.conj().T

    # cost function
    func = loc.cost_function

    # gradient function
    def grad(u):
        return loc.get_grad(u)

    # hessian diagonal function
    def hess_diag(u):
        return loc.gen_g_hop(u)[2]

    # hessian linear transformation function
    def hess_x(u, x):
        return loc.gen_g_hop(u)[1](x)

    # number of parameters
    n_param = (norb - 1) * norb // 2

    # call solver
    u = solver(unpack, func, grad, hess_diag, hess_x, n_param, "min", "cubic")
