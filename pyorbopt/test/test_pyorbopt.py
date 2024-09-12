from __future__ import annotations


import pytest
from pyscf import lo
from typing import TYPE_CHECKING

from pyorbopt.pyorbopt import solver

if TYPE_CHECKING:
    from pyscf import gto, scf


test_cases_solver = [
    ("h2o", lo.PM, "occ", "can", 4.101528606990511),
    ("h2o", lo.PM, "occ", "atomic", 4.101528606990511),
    ("h2o", lo.Boys, "occ", "can", 6.890557872084965),
    ("h2o", lo.Boys, "occ", "atomic", 6.890557872084965),
    ("h2o", lo.PM, "virt", "can", 7.101376595160159),
    ("h2o", lo.PM, "virt", "atomic", 7.101376595160159),
    ("h2o", lo.Boys, "virt", "can", 27.42565864888862),
    ("h2o", lo.Boys, "virt", "atomic", 27.42565864888862),
]


@pytest.mark.parametrize(
    argnames="system, orbs, space, guess, ref_cost_function",
    argvalues=test_cases_solver,
    ids=[
        "-".join([case[0], case[1].__name__, case[2], case[3]])
        for case in test_cases_solver
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
        norb = min(mol.nelec)
        mo_coeff = hf.mo_coeff[:, : min(mol.nelec)]
    elif space == "virt":
        norb = mol.nao - max(mol.nelec)
        mo_coeff = hf.mo_coeff[:, max(mol.nelec) :]

    # atomic initial guess
    if guess == "atomic":
        u = lo.boys.atomic_init_guess(mol, mo_coeff)
        mo_coeff = mo_coeff @ u

    # initialize pyscf localization object
    loc = orbs(mol, mo_coeff)

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

    # decide whether to maximize or minimize
    direction = "max" if isinstance(loc, lo.PM) else "min"

    # call solver
    u = solver(func, grad, hess_diag, hess_x, norb, direction)

    assert loc.cost_function(u) == pytest.approx(ref_cost_function)
