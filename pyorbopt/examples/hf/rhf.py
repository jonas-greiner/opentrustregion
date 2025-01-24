import numpy as np
from pyscf import gto, scf

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

# PySCF RHF calculation
hf = scf.RHF(mol).newton()
hf.conv_tol = 1.0e-10
hf.verbose = 4
hf.kernel()

# initialize object for second-order optimization
hf = scf.RHF(mol).newton()

# get one-electron Hamiltonian
h1e = hf._scf.get_hcore(mol)

# get one-electron overlap
s1e = hf._scf.get_ovlp(mol)

# get initial guess
dm = hf.get_init_guess(hf._scf.mol, hf.init_guess)

# generate hf potential
vhf = hf._scf.get_veff(mol, dm)

# get Fock matrix for initial guess
fock = hf.get_fock(h1e, s1e, vhf, dm)

# solve generalized eigenvalue problem
mo_energy, mo_coeff = hf.eig(fock, s1e)

# get orbital occupation
mo_occ = hf.get_occ(mo_energy, mo_coeff)

# get indices of all mixed occupation combinations
occidxa = mo_occ > 0
occidxb = mo_occ == 2
viridxa = ~occidxa
viridxb = ~occidxb
mask = (viridxa[:, None] & occidxa) | (viridxb[:, None] & occidxb)


# unpack matrix
def unpack(kappa):
    matrix = np.zeros(2 * (mol.nao,), dtype=np.float64)
    matrix[mask] = kappa
    return matrix - matrix.T


# energy function
def func(u):
    rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
    dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
    vhf = hf._scf.get_veff(mol, dm)
    return hf._scf.energy_tot(dm, h1e, vhf)


# energy and gradient function
def func_grad(u):
    rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
    dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
    vhf = hf._scf.get_veff(mol, dm)
    fock = hf.get_fock(h1e, s1e, vhf, dm)

    return hf._scf.energy_tot(dm, h1e, vhf), hf.get_grad(rot_mo_coeff, mo_occ, fock)


# energy, gradient, Hessian diagonal and Hessian linear transformation function
def func_grad_hdiag_hess_x(u):
    rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
    dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
    vhf = hf._scf.get_veff(mol, dm)
    fock = hf.get_fock(h1e, s1e, vhf, dm)

    g, hess_x, hdiag = hf.gen_g_hop(rot_mo_coeff, mo_occ, fock)

    return hf._scf.energy_tot(dm, h1e, vhf), g, hdiag, hess_x


# number of parameters
n_param = np.count_nonzero(mask)

# call solver
orbopt = OrbOpt(verbose=2)
u = orbopt.solver(
    unpack, func, func_grad, func_grad_hdiag_hess_x, n_param, "min", line_search="cubic"
)
