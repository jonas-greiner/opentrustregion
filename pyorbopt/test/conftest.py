from __future__ import annotations


import pytest
from pyscf import gto, scf
from typing import TYPE_CHECKING
from warnings import catch_warnings, simplefilter


if TYPE_CHECKING:
    from _pytest.fixtures import SubRequest


@pytest.fixture
def system(request: SubRequest) -> str:
    """
    this fixture stores the system string for other fixtures to access
    """
    return request.param


@pytest.fixture
def mol(system: str) -> gto.Mole:
    """
    this fixture constructs the mol object
    """
    if system == "h2o":
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

    return mol


@pytest.fixture
def hf(mol: gto.Mole) -> scf.RHF:
    """
    this fixture constructs the hf object and executes a hf calculation
    """
    hf = scf.RHF(mol)
    hf.conv_tol = 1.0e-10
    with catch_warnings():
        simplefilter("ignore")
        hf.kernel()

    return hf
