import pytest
import numpy as np
from ma1522 import Matrix
from ma1522.custom_types import (
    Shape,
    PartGen,
    ScalarFactor,
    PLU,
    RREF,
    VecDecomp,
    QR,
    PDP,
    SVD,
    NumSVD,
)


class TestCustomTypes:
    def test_shape_enum(self):
        assert Shape.SCALAR.value == "SCALAR"
        assert Shape.UPPER.value == "UPPER"
        assert Shape.LOWER.value == "LOWER"
        assert Shape.STRICT_UPPER.value == "STRICT_UPPER"
        assert Shape.STRICT_LOWER.value == "STRICT_LOWER"
        assert Shape.SYMMETRIC.value == "SYMMETRIC"

    def test_partgen(self):
        part_sol = Matrix([[1, 0], [0, 1]])
        gen_sol = Matrix([[0, 1], [1, 0]])
        pg = PartGen(part_sol, gen_sol)
        assert pg.part_sol == part_sol
        assert pg.gen_sol == gen_sol

    def test_scalarfactor(self):
        diag = Matrix([[1, 0], [0, 1]])
        full = Matrix([[1, 2], [3, 4]])
        sf = ScalarFactor(diag, full, "FD")
        assert sf.diag == diag
        assert sf.full == full
        assert sf.order == "FD"

    def test_plu(self):
        P = Matrix([[1, 0], [0, 1]])
        L = Matrix([[1, 0], [0, 1]])
        U = Matrix([[1, 0], [0, 1]])
        plu = PLU(P, L, U)
        assert plu.P == P
        assert plu.L == L
        assert plu.U == U

    def test_rref(self):
        rref_mat = Matrix([[1, 0], [0, 1]])
        pivots = (0, 1)
        rref = RREF(rref_mat, pivots)
        assert rref.rref == rref_mat
        assert rref.pivots == pivots

    def test_vecdecomp(self):
        proj = Matrix([[1], [0]])
        norm = Matrix([[0], [1]])
        vd = VecDecomp(proj, norm)
        assert vd.proj == proj
        assert vd.norm == norm

    def test_qr(self):
        Q = Matrix([[1, 0], [0, 1]])
        R = Matrix([[1, 0], [0, 1]])
        qr = QR(Q, R)
        assert qr.Q == Q
        assert qr.R == R

    def test_pdp(self):
        P = Matrix([[1, 0], [0, 1]])
        D = Matrix([[1, 0], [0, 1]])
        pdp = PDP(P, D)
        assert pdp.P == P
        assert pdp.D == D

    def test_svd(self):
        U = Matrix([[1, 0], [0, 1]])
        S = Matrix([[1, 0], [0, 1]])
        V = Matrix([[1, 0], [0, 1]])
        svd = SVD(U, S, V)
        assert svd.U == U
        assert svd.S == S
        assert svd.V == V

    def test_numsvd(self):
        U = np.array([[1, 0], [0, 1]])
        S = np.array([[1, 0], [0, 1]])
        V = np.array([[1, 0], [0, 1]])
        numsvd = NumSVD(U, S, V)
        assert numsvd.U.all() == U.all()
        assert numsvd.S.all() == S.all()
        assert numsvd.V.all() == V.all()
