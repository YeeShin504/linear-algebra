from symbolic import Matrix


def main():
    m = Matrix.create_rand_matrix(3, 3)
    return m


class Parent:
    def __init__(self) -> None:
        self.parent = "Parent"


class Child(Parent):
    def __init__(self) -> None:
        super().__init__()
        self.child = "Child"


if __name__ == "__main__":
    main()
