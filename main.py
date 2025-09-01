# Gaussian elimination implementation by Pedro Palacios for MAT-2233 at UT San Antonio
from colorama import Fore, init
from typing import Tuple, List

def main() -> None:
    matrix: List[List[float]] =  [[2.0, 3., 7.0, 3.0],
                            [8.0, 5.0, 3.0, 6.],
                            [1.0, 1.0, 4.0, 3.0],
                            [0.0, 3.0, 5.0, 8.0]]
    print("Matrix:")
    print_matrix(matrix)
    init(autoreset=True) # Initialize colorama for error messages and auto-resetting after each line

    non_zero_index_matrix = [(row_index, col_index) if col_value != 0 else (-1, -1) for row_index, row_value in enumerate(matrix) for col_index, col_value in enumerate(row_value)]
    nonzero_indicies: List[Tuple[int, int, bool]] = []
    pivot_coordinates: List[Tuple[int, int, bool]] = []

    non_zero_found = False
    first_pivot_coordinate = None
    # To check if there's at least one nonzero entry in each row excluding the last element, if the matrix is consistent, and saves the first potential pivot
    for i in range(len(matrix)):
        for (row, col) in non_zero_index_matrix:
            if row == i and col != len(matrix[i]): # If there exists a nonzero coordinate in the ith row that's not the last element
                non_zero_found = True
                equal_to_one = matrix[row][col] == 1

                # To add the coordinates of potential pivots, the first nonzero entry of the each row
                nonzero_indicies.append((row, col, equal_to_one))
                first_pivot_coordinate = nonzero_indicies[0]

                break # Implies there's at least 1 nonzero element excluding the last element
        if not non_zero_found and matrix[i][-1] != 0:
            # If there's only zeros in the coefficients and a nonzero in the row's last column—inconsistent
            print(Fore.RED + "ERROR: No nonzero coefficient found, yet there's still a nonzero constant. This matrix is inconsistent.")
            return
    print(Fore.GREEN + "SUCCESS: The augmented matrix is consistent. Proceeding with creating pivots and performing row operations for RREF.")


    # To populate the pivot coordinates list with the expected coordinates of each pivot
    if first_pivot_coordinate is not None:
        new_pivot_row, new_pivot_column, new_pivot_equal_to_one = first_pivot_coordinate
        for i in range(first_pivot_coordinate[0], len(matrix)): # Start at the row of the first pivot and end at the last row of matrix
            if matrix[new_pivot_row][new_pivot_column] == 0 and new_pivot_column != (len(matrix[new_pivot_row]) - 2):
                new_pivot_column += 1
            elif i != (len(matrix) - 1) and new_pivot_column != (len(matrix[new_pivot_row]) - 1):
                    pivot_coordinates.append((new_pivot_row, new_pivot_column, new_pivot_equal_to_one))
                    new_pivot_row, new_pivot_column = new_pivot_row + 1, new_pivot_column + 1
                    new_pivot_equal_to_one = True if matrix[new_pivot_row][new_pivot_column] == 1 else False


    else:
        print(Fore.YELLOW + "WARNING: No pivot found in the matrix.")

    # To scale at expected pivot coordinates and eliminate nonzeros above and below pivots through row replacement
    while len(pivot_coordinates) > 0:
        for (pivot_row, pivot_column, pivot_equal_to_one) in pivot_coordinates:
            # To create a pivot equal to 1 at each coordinate such that on the following row to the right, the element is also equal to 1 for each row that's not the last one
            for i in range(pivot_row, len(matrix)):
                pivot_value = matrix[pivot_row][pivot_column]
                row_scaling(matrix, i, 1 / float(pivot_value), float(pivot_value))  # To change the pivot to 1

            print(f"Matrix after scaling:")
            print_matrix(matrix)

            # Search for nonzeros below each pivot and eliminate with pivot
            if pivot_row != len(matrix) - 1: # Verify the pivot is not in the last row of the matrix
                for i in range(pivot_row + 1, len(matrix)):
                    # Find any nonzeros below the pivot to eliminate through row replacement
                    value = float(matrix[i][pivot_column])
                    if value != 0:
                       row_replacement(matrix, i, -value, pivot_row)
                    else:
                       continue
                if pivot_row > 0:
                    for i in range(pivot_row - 1, -1, -1):
                        # Find any nonzeros above the pivot to eliminate through row replacement
                        value = float(matrix[i][pivot_column])
                        if value != 0:
                            row_replacement(matrix, i, -value, pivot_row)
                        else:
                            continue

            pivot_coordinates.remove((pivot_row, pivot_column, pivot_equal_to_one))
            print("Matrix after row replacement:")
            print_matrix(matrix)
    print(Fore.GREEN + "\nSUCCESS: Transformed augmented matrix into reduced row echelon form (RREF).")
    print_matrix(matrix)

def row_scaling(matrix: List[List[float]], row_index: int, num: float, pivot_value: float = None) -> None: # The pivot is optional and only for printing the multiplier
    for index, value in enumerate(matrix[row_index]):
        matrix[row_index][index] = value * num
    if pivot_value is not None:
        print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by 1/{pivot_value}.")
    else:
        print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by {num}")

def row_replacement(matrix: List[List[float]], row_index: int, multiplier: float, pivot_row: int) -> None:
    if multiplier == 0:
        print(Fore.RED + f"ERROR: Cannot perform row replacement with a multiplier equal to 0!")
    for index, value in enumerate(matrix[row_index]):
        matrix[row_index][index] = value + multiplier * matrix[pivot_row][index] # Current value + multiplier * value above in pivot row

    print(Fore.GREEN + f"SUCCESS: Replaced row {row_index} with R{row_index} + {multiplier}*R{pivot_row}")

def row_exchange(matrix: List[List[float]]) -> None:
    pass

def print_matrix(matrix: List[List[float]]):
    for row in matrix:
        print(row)

if __name__ == "__main__":
    main()