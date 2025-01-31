import numpy as np
from pyscf import gto, scf, lo, soscf

from pyorbopt.pyorbopt import OrbOpt

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
    def unpack(kappa: np.ndarray) -> np.ndarray:
        matrix = np.zeros(2 * (norb,), dtype=np.float64)
        idx = np.tril_indices(norb, -1)
        matrix[idx] = kappa
        return matrix - matrix.conj().T

    # cost function
    def func(kappa: np.ndarray) -> float:
        u = soscf.ciah.expmat(unpack(kappa))
        return loc.cost_function(u)

    # cost function, gradient, Hessian diagonal and Hessian linear transformation
    # function
    def update_orbs(kappa):
        u = soscf.ciah.expmat(unpack(kappa))
        func = loc.cost_function(u)
        grad, hess_x, hdiag = loc.gen_g_hop(u)
        loc.mo_coeff = loc.mo_coeff @ u
        return func, grad, hdiag, hess_x

    # number of parameters
    n_param = (norb - 1) * norb // 2

    # call solver
    orbopt = OrbOpt(verbose=2)
    u = orbopt.solver(func, update_orbs, n_param, line_search=True)
