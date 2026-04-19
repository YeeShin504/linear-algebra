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
    RREFCase,
)


class TestCustomTypes:
    def test_shape_enum(self):
        assert Shape.SCALAR.value == "SCALAR"
        assert Shape.UPPER.value == "UPPER"
        assert Shape.LOWER.value == "LOWER"
        assert Shape.STRICT_UPPER.value == "STRICT_UPPER"
        assert Shape.STRICT_LOWER.value == "STRICT_LOWER"
        assert Shape.SYMMETRIC.value == "SYMMETRIC"

    def test_printable_base(self):
        mat = Matrix([[1, 2], [3, 4]])
        pivots = (0, 1)
        rref = RREF(mat, pivots)
        
        # Iteration
        fields = list(rref)
        assert fields[0] == mat
        assert fields[1] == pivots
        
        # Get/Set
        assert rref[0] == mat
        new_mat = Matrix.eye(2)
        rref[0] = new_mat
        assert rref.rref == new_mat
        
        # Eval / Evalf
        assert rref.eval() == new_mat
        assert isinstance(rref.evalf(), Matrix)

    def test_partgen_methods(self):
        part = Matrix([[1], [0]])
        gen = Matrix([[0], [1]])
        pg = PartGen(part, gen)
        assert pg.eval() == Matrix([[1], [1]])
        assert "\\left(" in pg._latex()

    def test_scalarfactor_methods(self):
        diag = Matrix([[2, 0], [0, 2]])
        full = Matrix([[1, 2], [3, 4]])
        sf = ScalarFactor(diag, full, "DF")
        assert sf.eval() == Matrix([[2, 4], [6, 8]])
        assert sf._latex().startswith(diag._latex())
        
        sf_fd = ScalarFactor(diag, full, "FD")
        assert sf_fd.eval() == Matrix([[2, 4], [6, 8]])
        assert sf_fd._latex().endswith(diag._latex())

    def test_plu_methods(self):
        P, L, U = Matrix.eye(2), Matrix([[1, 0], [2, 1]]), Matrix([[3, 4], [0, 5]])
        plu = PLU(P, L, U)
        assert plu.eval() == L @ U
        assert plu._latex().count("array") == 6

    def test_vecdecomp_methods(self):
        proj, norm = Matrix([[1], [0]]), Matrix([[0], [1]])
        vd = VecDecomp(proj, norm)
        assert vd.eval() == Matrix([[1], [1]])

    def test_qr_methods(self):
        Q, R = Matrix.eye(2), Matrix([[1, 2], [0, 3]])
        qr = QR(Q, R)
        assert qr.eval() == R
        assert qr._latex().count("array") == 4

    def test_pdp_methods(self):
        # Regular case
        P = Matrix([[1, 1], [0, 1]])
        D = Matrix([[2, 0], [0, 3]])
        pdp = PDP(P, D)
        assert pdp.eval() == P @ D @ P.inv()
        assert pdp._latex().count("array") == 6
        
        # Singular case (exception fallback)
        P_sing = Matrix([[1, 1], [1, 1]])
        pdp_sing = PDP(P_sing, D)
        assert "P inverse does not exist" in pdp_sing._latex()

    def test_svd_methods(self):
        U, S, V = Matrix.eye(2), Matrix([[2, 0], [0, 1]]), Matrix.eye(2)
        svd = SVD(U, S, V)
        assert svd.eval() == S
        assert svd._latex().count("array") == 6

    def test_numsvd_methods(self):
        U = np.eye(2)
        S = np.diag([2, 1])
        V = np.eye(2)
        nsvd = NumSVD(U, S, V)
        assert np.allclose(nsvd.eval(), S)
        assert "NumSVD" in repr(nsvd)

    def test_rref_case_methods(self):
        rc = RREFCase(
            conditions={"a": 0},
            excluded=[{"b": 0}],
            rref=Matrix.eye(2),
            pivots=(0, 1),
            free_params=0,
            is_consistent=True
        )
        assert rc.eval() == Matrix.eye(2)
        assert "RREFCase" in rc._latex()

    def test_repr_latex(self):
        rref = RREF(Matrix.eye(2), (0, 1))
        assert rref._repr_latex_().startswith("$")
        assert rref._repr_latex_().endswith("$")
