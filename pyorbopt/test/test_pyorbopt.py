from __future__ import annotations


import pytest
import numpy as np
import scipy as sc
from pyscf import scf, lo
from typing import TYPE_CHECKING

from pyorbopt.pyorbopt import OrbOpt

if TYPE_CHECKING:
    from typing import Tuple, Callable

    from pyscf import gto, scf


test_cases_solver_loc = [
    ("h2o", lo.PM, "occ", "can", 4.101528606990511),
    ("h2o", lo.PM, "occ", "atomic", 4.101528606990511),
    ("h2o", lo.Boys, "occ", "can", 6.890557872084965),
    ("h2o", lo.Boys, "occ", "atomic", 6.890557872084965),
    ("h2o", lo.PM, "virt", "can", 7.101376595160159),
    ("h2o", lo.PM, "virt", "atomic", 7.101376595160159),
    ("h2o", lo.Boys, "virt", "can", 27.42565864888862),
    ("h2o", lo.Boys, "virt", "atomic", 27.42565864888862),
]

test_cases_solver_hf = [
    ("h2o", scf.hf_symm.RHF, -75.98399811804387),
    ("h2o+", scf.hf_symm.ROHF, -75.57838419175835),
    ("h2o+", scf.uhf_symm.UHF, -75.5805066388706),
]


@pytest.mark.parametrize(
    argnames="system, orbs, space, guess, ref_cost_function",
    argvalues=test_cases_solver_loc,
    ids=[
        "-".join([case[0], case[1].__name__, case[2], case[3]])
        for case in test_cases_solver_loc
    ],
    indirect=["system"],
)
def test_solver_loc(
    mol: gto.Mole,
    hf: scf.RHF,
    orbs: lo.OrbitalLocalizer,
    space: str,
    guess: str,
    ref_cost_function: float,
) -> None:
    """
    this function tests the solver for orbital localization
    """
    # canonical initial guess
    if space == "occ":
        mo_coeff = hf.mo_coeff[:, : min(mol.nelec)]
    elif space == "virt":
        mo_coeff = hf.mo_coeff[:, max(mol.nelec) :]

    # number of orbitals
    norb = mo_coeff.shape[1]

    # atomic initial guess
    if guess == "atomic":
        u = lo.boys.atomic_init_guess(mol, mo_coeff)
        mo_coeff = mo_coeff @ u

    # initialize pyscf localization object
    loc = orbs(mol, mo_coeff)

    # unpack matrix
    def unpack(kappa: np.ndarray) -> np.ndarray:
        matrix = np.zeros(2 * (norb,), dtype=np.float64)
        idx = np.tril_indices(norb, -1)
        matrix[idx] = kappa
        return matrix - matrix.conj().T

    # cost function
    def func(kappa: np.ndarray) -> float:
        u = sc.linalg.expm(unpack(kappa))
        if isinstance(loc, lo.PM):
            return -loc.cost_function(u)
        else:
            return loc.cost_function(u)

    # cost function, gradient, Hessian diagonal and Hessian linear transformation
    # function
    def update_orbs(
        kappa: np.ndarray,
    ) -> Tuple[float, np.ndarray, np.ndarray, Callable[[np.ndarray], np.ndarray]]:
        u = sc.linalg.expm(unpack(kappa))
        func = loc.cost_function(u)
        grad, hess_x, hdiag = loc.gen_g_hop(u)
        loc.mo_coeff = mo_coeff @ u
        if isinstance(loc, lo.PM):
            return -func, -grad, -hdiag, lambda x: -hess_x(x)
        else:
            return func, grad, hdiag, hess_x

    # number of parameters
    n_param = (norb - 1) * norb // 2

    # call solver
    orbopt = OrbOpt(verbose=3)
    u = orbopt.solver(func, update_orbs, n_param)

    assert loc.cost_function(u) == pytest.approx(ref_cost_function)


@pytest.mark.parametrize(
    argnames="system, scf_class, ref_energy",
    argvalues=test_cases_solver_hf,
    ids=["-".join([case[0], case[1].__name__]) for case in test_cases_solver_hf],
    indirect=["system"],
)
def test_solver_hf(mol: gto.Mole, scf_class: scf.SCF, ref_energy: float) -> None:
    """
    this function tests the solver for orbital localization
    """
    hf = scf_class(mol).newton()
    hf.conv_tol = 1.0e-10

    hf.verbose = 4
    hf.kernel()

    class HFWrapper:
        def __init__(self, mol):
            # initialize mol object
            self.mol = mol

            # initialize object for second-order optimization
            self.hf = scf_class(mol).newton()

            # get one-electron Hamiltonian
            self.h1e = hf._scf.get_hcore(mol)

            # get one-electron overlap
            self.s1e = hf._scf.get_ovlp(mol)

            # get initial guess
            dm = hf.get_init_guess(hf._scf.mol, hf.init_guess)

            # generate hf potential
            vhf = hf._scf.get_veff(mol, dm)

            # get Fock matrix for initial guess
            fock = hf.get_fock(h1e=self.h1e, s1e=self.s1e, vhf=vhf, dm=dm)

            # solve generalized eigenvalue problem
            mo_energy, self.mo_coeff = hf.eig(fock, self.s1e)

            # get orbital occupation
            self.mo_occ = hf.get_occ(mo_energy, self.mo_coeff)

            # get indices of all mixed occupation combinations
            if self.mo_occ.ndim == 1:
                occidxa = self.mo_occ > 0
                occidxb = self.mo_occ == 2
                viridxa = ~occidxa
                viridxb = ~occidxb
                self.mask = (viridxa[:, None] & occidxa) | (viridxb[:, None] & occidxb)
            else:
                occidxa = self.mo_occ[0] == 1
                occidxb = self.mo_occ[1] == 1
                viridxa = ~occidxa
                viridxb = ~occidxb
                self.maska = viridxa[:, None] & occidxa
                self.maskb = viridxb[:, None] & occidxb

        # unpack matrix
        def unpack(self, kappa: np.ndarray) -> np.ndarray:
            if self.mo_occ.ndim == 1:
                matrix = np.zeros(2 * (self.mol.nao,), dtype=np.float64)
                matrix[self.mask] = kappa
            else:
                matrix = np.zeros(2 * (2 * self.mol.nao,), dtype=np.float64)
                matrix[: self.mol.nao, : mol.nao][self.maska] = kappa[
                    : np.count_nonzero(self.maska)
                ]
                matrix[self.mol.nao :, self.mol.nao :][self.maskb] = kappa[
                    np.count_nonzero(self.maska) :
                ]
            return matrix - matrix.T

        # energy function
        def func(self, kappa: np.ndarray) -> float:
            u = sc.linalg.expm(self.unpack(kappa))
            if self.mo_occ.ndim == 1:
                rot_mo_coeff = self.hf.rotate_mo(self.mo_coeff, u)
            else:
                rot_mo_coeff = self.hf.rotate_mo(
                    self.mo_coeff,
                    (
                        u[: self.mol.nao, : self.mol.nao],
                        u[self.mol.nao :, self.mol.nao :],
                    ),
                )
            dm = self.hf.make_rdm1(rot_mo_coeff, self.mo_occ)
            vhf = self.hf._scf.get_veff(mol, dm=dm)
            return self.hf._scf.energy_tot(dm=dm, h1e=self.h1e, vhf=vhf)

        # energy, gradient, Hessian diagonal and Hessian linear transformation function
        def update_orbs(
            self,
            kappa: np.ndarray,
        ) -> Tuple[float, np.ndarray, np.ndarray, Callable[[np.ndarray], np.ndarray]]:
            u = sc.linalg.expm(self.unpack(kappa))
            if self.mo_occ.ndim == 1:
                self.mo_coeff = self.hf.rotate_mo(self.mo_coeff, u)
            else:
                self.mo_coeff = self.hf.rotate_mo(
                    self.mo_coeff,
                    (
                        u[: self.mol.nao, : self.mol.nao],
                        u[self.mol.nao :, self.mol.nao :],
                    ),
                )
            dm = self.hf.make_rdm1(self.mo_coeff, self.mo_occ)
            vhf = self.hf._scf.get_veff(mol, dm=dm)
            fock = self.hf.get_fock(h1e=self.h1e, s1e=self.s1e, vhf=vhf, dm=dm)
            g, hess_x, hdiag = self.hf.gen_g_hop(self.mo_coeff, self.mo_occ, fock)
            return (
                self.hf._scf.energy_tot(dm=dm, h1e=self.h1e, vhf=vhf),
                g,
                hdiag,
                hess_x,
            )

    hf_wrapper = HFWrapper(mol)

    # number of parameters
    if hf_wrapper.mo_occ.ndim == 1:
        n_param = np.count_nonzero(hf_wrapper.mask)
    else:
        n_param = np.count_nonzero(hf_wrapper.maska) + np.count_nonzero(
            hf_wrapper.maskb
        )

    # call solver
    orbopt = OrbOpt(verbose=3)
    orbopt.solver(
        hf_wrapper.func,
        hf_wrapper.update_orbs,
        n_param,
        conv_tol=1e-8,
        n_random_trial_vectors=2,
    )
    print(hf_wrapper.func(np.zeros(n_param)))
    print(hf_wrapper.update_orbs(np.zeros(n_param))[0])
    print(np.linalg.norm(hf_wrapper.update_orbs(np.zeros(n_param))[1]))

    assert hf_wrapper.func(np.zeros(n_param)) == pytest.approx(ref_energy)
