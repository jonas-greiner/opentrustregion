from __future__ import annotations


import pytest
import numpy as np
from pyscf import scf, lo
from typing import TYPE_CHECKING

from pyorbopt.pyorbopt import solver

if TYPE_CHECKING:
    from pyscf import gto, scf
    from typing import Callable


test_cases_solver_loc = [
    ("h2o", lo.PM, "occ", "can", "golden", 4.101528606990511),
    ("h2o", lo.PM, "occ", "can", "brent", 4.101528606990511),
    ("h2o", lo.PM, "occ", "can", "cubic", 4.101528606990511),
    ("h2o", lo.PM, "occ", "atomic", "golden", 4.101528606990511),
    ("h2o", lo.PM, "occ", "atomic", "brent", 4.101528606990511),
    ("h2o", lo.PM, "occ", "atomic", "cubic", 4.101528606990511),
    ("h2o", lo.Boys, "occ", "can", "golden", 6.890557872084965),
    ("h2o", lo.Boys, "occ", "can", "brent", 6.890557872084965),
    ("h2o", lo.Boys, "occ", "can", "cubic", 6.890557872084965),
    ("h2o", lo.Boys, "occ", "atomic", "golden", 6.890557872084965),
    ("h2o", lo.Boys, "occ", "atomic", "brent", 6.890557872084965),
    ("h2o", lo.Boys, "occ", "atomic", "cubic", 6.890557872084965),
    ("h2o", lo.PM, "virt", "can", "golden", 7.101376595160159),
    ("h2o", lo.PM, "virt", "can", "brent", 7.101376595160159),
    ("h2o", lo.PM, "virt", "can", "cubic", 7.101376595160159),
    ("h2o", lo.PM, "virt", "atomic", "golden", 7.101376595160159),
    ("h2o", lo.PM, "virt", "atomic", "brent", 7.101376595160159),
    ("h2o", lo.PM, "virt", "atomic", "cubic", 7.101376595160159),
    ("h2o", lo.Boys, "virt", "can", "golden", 27.42565864888862),
    ("h2o", lo.Boys, "virt", "can", "brent", 27.42565864888862),
    ("h2o", lo.Boys, "virt", "can", "cubic", 27.42565864888862),
    ("h2o", lo.Boys, "virt", "atomic", "golden", 27.42565864888862),
    ("h2o", lo.Boys, "virt", "atomic", "brent", 27.42565864888862),
    ("h2o", lo.Boys, "virt", "atomic", "cubic", 27.42565864888862),
]

test_cases_solver_hf = [
    ("h2o", scf.hf_symm.RHF, "golden", -75.98399811804387),
    ("h2o", scf.hf_symm.RHF, "brent", -75.98399811804387),
    ("h2o", scf.hf_symm.RHF, "cubic", -75.98399811804387),
    ("h2o+", scf.hf_symm.ROHF, "golden", -75.57838419175835),
    ("h2o+", scf.hf_symm.ROHF, "brent", -75.57838419175835),
    ("h2o+", scf.hf_symm.ROHF, "cubic", -75.57838419175835),
    ("h2o+", scf.uhf_symm.UHF, "golden", -75.5805066388706),
    ("h2o+", scf.uhf_symm.UHF, "brent", -75.5805066388706),
    ("h2o+", scf.uhf_symm.UHF, "cubic", -75.5805066388706),
]


@pytest.mark.parametrize(
    argnames="system, orbs, space, guess, line_search, ref_cost_function",
    argvalues=test_cases_solver_loc,
    ids=[
        "-".join([case[0], case[1].__name__, case[2], case[3], case[4]])
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
    line_search: str,
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
    func = loc.cost_function

    # gradient function
    def grad(u: np.ndarray) -> np.ndarray:
        return loc.get_grad(u)

    # hessian diagonal function
    def hess_diag(u: np.ndarray) -> np.ndarray:
        return loc.gen_g_hop(u)[2]

    # hessian linear transformation function
    def hess_x(u: np.ndarray, x: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
        return loc.gen_g_hop(u)[1](x)

    # number of parameters
    n_param = (norb - 1) * norb // 2

    # decide whether to maximize or minimize
    direction = "max" if isinstance(loc, lo.PM) else "min"

    # call solver
    u = solver(unpack, func, grad, hess_diag, hess_x, n_param, direction, line_search)

    assert loc.cost_function(u) == pytest.approx(ref_cost_function)


@pytest.mark.parametrize(
    argnames="system, scf_class, line_search, ref_energy",
    argvalues=test_cases_solver_hf,
    ids=[
        "-".join([case[0], case[1].__name__, case[2]]) for case in test_cases_solver_hf
    ],
    indirect=["system"],
)
def test_solver_hf(
    mol: gto.Mole, scf_class: scf.SCF, line_search: str, ref_energy: float
) -> None:
    """
    this function tests the solver for orbital localization
    """
    hf = scf_class(mol).newton()
    hf.conv_tol = 1.0e-10

    hf.verbose = 4
    hf.kernel()

    # initialize object for second-order optimization
    hf = scf_class(mol).newton()

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
    if mo_occ.ndim == 1:
        occidxa = mo_occ > 0
        occidxb = mo_occ == 2
        viridxa = ~occidxa
        viridxb = ~occidxb
        mask = (viridxa[:, None] & occidxa) | (viridxb[:, None] & occidxb)
    else:
        occidxa = mo_occ[0] == 1
        occidxb = mo_occ[1] == 1
        viridxa = ~occidxa
        viridxb = ~occidxb
        maska = viridxa[:, None] & occidxa
        maskb = viridxb[:, None] & occidxb

    # unpack matrix
    def unpack(kappa: np.ndarray) -> np.ndarray:
        if mo_occ.ndim == 1:
            matrix = np.zeros(2 * (mol.nao,), dtype=np.float64)
            matrix[mask] = kappa
        else:
            matrix = np.zeros(2 * (2 * mol.nao,), dtype=np.float64)
            matrix[: mol.nao, : mol.nao][maska] = kappa[: np.count_nonzero(maska)]
            matrix[mol.nao :, mol.nao :][maskb] = kappa[np.count_nonzero(maska) :]
        return matrix - matrix.T

    # energy function
    def func(u: np.ndarray) -> np.ndarray:
        if mo_occ.ndim == 2:
            u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
        rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
        dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
        vhf = hf._scf.get_veff(mol, dm)
        return hf._scf.energy_tot(dm, h1e, vhf)

    # gradient function
    def grad(u: np.ndarray) -> np.ndarray:
        if mo_occ.ndim == 2:
            u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
        rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
        dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
        vhf = hf._scf.get_veff(mol, dm)
        fock_ao = hf.get_fock(h1e, s1e, vhf, dm)
        return hf.get_grad(rot_mo_coeff, mo_occ, fock_ao)

    # hessian diagonal function
    def hess_diag(u: np.ndarray) -> np.ndarray:
        if mo_occ.ndim == 2:
            u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
        rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
        dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
        vhf = hf._scf.get_veff(mol, dm)
        fock_ao = hf.get_fock(h1e, s1e, vhf, dm)
        return hf.gen_g_hop(rot_mo_coeff, mo_occ, fock_ao)[2]

    # hessian linear transformation function
    def hess_x(u: np.ndarray, x: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
        if mo_occ.ndim == 2:
            u = (u[: mol.nao, : mol.nao], u[mol.nao :, mol.nao :])
        rot_mo_coeff = hf.rotate_mo(mo_coeff, u)
        dm = hf.make_rdm1(rot_mo_coeff, mo_occ)
        vhf = hf._scf.get_veff(mol, dm)
        fock_ao = hf.get_fock(h1e, s1e, vhf, dm)
        return hf.gen_g_hop(rot_mo_coeff, mo_occ, fock_ao)[1](x)

    # number of parameters
    if mo_occ.ndim == 1:
        n_param = np.count_nonzero(mask)
    else:
        n_param = np.count_nonzero(maska) + np.count_nonzero(maskb)

    # call solver
    u = solver(unpack, func, grad, hess_diag, hess_x, n_param, "min", line_search)

    assert func(u) == pytest.approx(ref_energy)
