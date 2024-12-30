from multiprocessing import Pool

class Shell:
    def __init__(self) -> None:
        self.y0 = [1, 2, 3, 4]
        self.y_NK_Cell = [1, 2, 3, 4]
        self.tD = 4.0

    def Sol_Lysis(self, tD, y0, y_NK_Cell):
        return sum(y0) + sum(y_NK_Cell) + tD

    @staticmethod
    def worker(args):
        # Unpack the arguments
        self, E = args
        print(self.y0)
        # Multiply y_NK_Cell with E
        y_NK_Cell_multiplied = [y * E for y in self.y_NK_Cell]
        # Call the Sol_Lysis method
        yF = self.Sol_Lysis(self.tD, self.y0, y_NK_Cell_multiplied)
        return yF

    def test(self):
        Effector = [100000, 50000, 25000, 12500, 6200]
        # Create a tuple of arguments for each E in Effector
        args_list = [(self, E) for E in Effector]
        with Pool() as pool:
            Specific_lysis = pool.map(Shell.worker, args_list)
        print(Specific_lysis)

if __name__ == "__main__":
    a = Shell()
    a.test()