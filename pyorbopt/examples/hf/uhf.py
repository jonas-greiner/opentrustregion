import numpy as np
from pyscf import gto, scf

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
mol.charge = 1
mol.spin = 1

# PySCF UHF calculation
hf = scf.UHF(mol).newton()
hf.conv_tol = 1.0e-10
hf.verbose = 4
hf.kernel()

# initialize object for second-order optimization
hf = scf.UHF(mol).newton()

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
occidxa = mo_occ[0] == 1
occidxb = mo_occ[1] == 1
viridxa = ~occidxa
viridxb = ~occidxb
maska = viridxa[:, None] & occidxa
maskb = viridxb[:, None] & occidxb


# unpack matrix
def unpack(kappa):
    matrix = np.zeros(2 * (2 * mol.nao,), dtype=np.float64)
    matrix[: mol.nao, : mol.nao][maska] = kappa[: np.count_nonzero(maska)]
    matrix[mol.nao :, mol.nao :][maskb] = kappa[np.count_nonzero(maska) :]
    return matrix - matrix.T


# energy function
def func(u):
    if mo_occ.ndim == 2:
        u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
    rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
    dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
    vhf = hf._scf.get_veff(mol, dm)
    return hf._scf.energy_tot(dm, h1e, vhf)


# gradient function
def grad(u):
    if mo_occ.ndim == 2:
        u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
    rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
    dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
    vhf = hf._scf.get_veff(mol, dm)
    fock_ao = hf.get_fock(h1e, s1e, vhf, dm)
    return hf.get_grad(rot_mo_coeff, mo_occ, fock_ao)


# hessian diagonal function
def hess_diag(u):
    if mo_occ.ndim == 2:
        u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
    rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
    dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
    vhf = hf._scf.get_veff(mol, dm)
    fock_ao = hf.get_fock(h1e, s1e, vhf, dm)
    return hf.gen_g_hop(rot_mo_coeff, mo_occ, fock_ao)[2]


# hessian linear transformation function
def hess_x(u, x):
    if mo_occ.ndim == 2:
        u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
    rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
    dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
    vhf = hf._scf.get_veff(mol, dm)
    fock_ao = hf.get_fock(h1e, s1e, vhf, dm)
    return hf.gen_g_hop(rot_mo_coeff, mo_occ, fock_ao)[1](x)


# number of parameters
n_param = np.count_nonzero(maska) + np.count_nonzero(maskb)

# call solver
u = solver(unpack, func, grad, hess_diag, hess_x, n_param, "min", "cubic")
