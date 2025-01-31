import numpy as np
from pyscf import gto, scf, soscf

from pyorbopt.pyorbopt import OrbOpt

mol = gto.Mole()
mol.build(
    atom="""
    O	 0.0000000	 0.0000000	 0.0000000
    H	 0.7569685	 0.0000000	-0.5858752
    H	-0.7569685	 0.0000000	-0.5858752
""",
    basis="631g",
    symmetry=False,
    verbose=0,
)

# PySCF RHF calculation
hf = scf.RHF(mol).newton()
hf.conv_tol = 1.0e-10
hf.verbose = 4
hf.kernel()


class HFWrapper:
    def __init__(self, mol):
        # initialize mol object
        self.mol = mol

        # initialize object for second-order optimization
        self.hf = scf.RHF(mol).newton()

        # get one-electron Hamiltonian
        self.h1e = self.hf._scf.get_hcore(mol)

        # get one-electron overlap
        self.s1e = self.hf._scf.get_ovlp(mol)

        # get initial guess
        dm = self.hf.get_init_guess(self.hf._scf.mol, self.hf.init_guess)

        # generate hf potential
        vhf = self.hf._scf.get_veff(mol, dm)

        # get Fock matrix for initial guess
        fock = self.hf.get_fock(self.h1e, self.s1e, vhf, dm)

        # solve generalized eigenvalue problem
        mo_energy, self.mo_coeff = self.hf.eig(fock, self.s1e)

        # get orbital occupation
        self.mo_occ = self.hf.get_occ(mo_energy, self.mo_coeff)

        # get indices of all mixed occupation combinations
        occidxa = self.mo_occ > 0
        occidxb = self.mo_occ == 2
        viridxa = ~occidxa
        viridxb = ~occidxb
        self.mask = (viridxa[:, None] & occidxa) | (viridxb[:, None] & occidxb)

    # unpack matrix
    def unpack(self, kappa):
        matrix = np.zeros(2 * (self.mol.nao,), dtype=np.float64)
        matrix[self.mask] = kappa
        return matrix - matrix.T

    # energy function
    def func(self, kappa):
        u = soscf.ciah.expmat(self.unpack(kappa))
        rot_mo_coeff = self.hf.rotate_mo(self.mo_coeff, u)
        dm = self.hf.make_rdm1(rot_mo_coeff, self.mo_occ)
        vhf = self.hf._scf.get_veff(mol, dm)
        return self.hf._scf.energy_tot(dm, self.h1e, vhf)

    # energy, gradient, Hessian diagonal and Hessian linear transformation function
    def update_orbs(self, kappa):
        u = soscf.ciah.expmat(self.unpack(kappa))
        self.mo_coeff = self.hf.rotate_mo(self.mo_coeff, u)
        dm = self.hf.make_rdm1(self.mo_coeff, self.mo_occ)
        vhf = self.hf._scf.get_veff(mol, dm)
        fock = self.hf.get_fock(self.h1e, self.s1e, vhf, dm)
        g, hess_x, hdiag = self.hf.gen_g_hop(self.mo_coeff, self.mo_occ, fock)
        return (self.hf._scf.energy_tot(dm, self.h1e, vhf), g, hdiag, hess_x)


hf_wrapper = HFWrapper(mol)

# number of parameters
n_param = np.count_nonzero(hf_wrapper.mask)

# call solver
orbopt = OrbOpt(verbose=2)
u = orbopt.solver(hf_wrapper.func, hf_wrapper.update_orbs, n_param)
