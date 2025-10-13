from gaussianelimination import GaussianElimination
from manim import *
from typing import List

class GaussianVisualizations(Scene):
    def construct(self):
        Matrix = List[List[float]]
        # matrix: Matrix = [[0, 1, 5, 3],
        #                             [1, 4, 4, -3],
        #                             [2, 6, 3, -2]]

        matrix: Matrix = [[1, 4, -5, 1, 2],
                  [2, 5, -4, -1, 4],
                  [-3, -9, 9, 2, 2],
                  [1, 7, 3, 6, 9]]

        # matrix: Matrix = [[2, 1, 4],
        #           [2, -1, 0]]

        # matrix: Matrix = [[0, 0, 6, 2],
        #           [0, 1, 4, -7],
        #           [1, 1, 10, -5]]

        # matrix: Matrix = [[0, 0, -5, 0, -6, 3],
        #           [0, 2, 8, -1, 0, 2],
        #           [0, 0, 0, 0, 2, 0],
        #           [0, 0, 0, 0, 0, 0]]

        # matrix: Matrix = [[1, 0, -6, 0, 3],
        #           [3, 4, 8, -1, 2],
        #           [0, 9, -4, 2, -1],
        #           [0, 7, -7, 18, 0],
        #           [1, 2, -5, 0, 3],
        #           [0, 1, 8, -1, 2],
        #           [8, 9, 13, 0, 0]]

        starting_decimal_matrix: DecimalMatrix = DecimalMatrix(matrix)

        # Run gaussian elimination on matrix
        gaussianelimination = GaussianElimination(matrix)
        matrix: Matrix = gaussianelimination.get_matrix()
        text = Text("Gaussian Elimination", font="Times New Roman").move_to((0, 2.5, 0))
        decimal_matrix: DecimalMatrix = DecimalMatrix(matrix)

        # Run animations on starting and final matrix
        self.play(FadeIn(starting_decimal_matrix))
        self.play(FadeOut(starting_decimal_matrix, shift=(-3, 0, 0)))
        group = Group(starting_decimal_matrix, decimal_matrix)
        group.arrange_in_grid(buff=2).scale(0.85)
        arrow = Arrow(start=LEFT, end=RIGHT).next_to(starting_decimal_matrix, buff=0.12)
        self.play(FadeIn(starting_decimal_matrix), FadeIn(arrow), FadeIn(text),FadeIn(decimal_matrix))
