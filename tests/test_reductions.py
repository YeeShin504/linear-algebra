import pytest


from symbolic import Matrix


class TestAugLine:
    def test1(self):
        mat = Matrix([[1, 2], [3, 4]], aug_pos=0)
        mat2 = mat.aug_line(1)
        assert hasattr(mat2, "_aug_pos")
        assert mat2._aug_pos == set([1, 0])
