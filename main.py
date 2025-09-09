# Gaussian elimination implementation by Pedro Palacios for MAT-2233 at UT San Antonio
from colorama import Fore, init
from typing import Tuple, List

def main() -> None:
    global pivot_coordinates
    # matrix: List[List[float]] = [[0, 1, 5, 3],
    #                             [1, 4, 4, -3],
    #                             [2, 6, 3, -2]]

    matrix = [[1, 4, -5, 1, 2],
              [2, 5, -4, -1, 4],
              [-3, -9, 9, 2, 2]]

    # matrix = [[2, 1, 4],
    #           [2, -1, 0]]

    # matrix = [[0, 0, 6, 2],
    #           [0, 1, 4, -7],
    #           [1, 1, 10, -5]]

    # matrix = [[0, 0, -5, 0, -6, 3],
    #           [0, 2, 8, -1, 0, 2],
    #           [0, 0, 0, 0, 2, 0],
    #           [0, 0, 0, 0, 0, 0]]

    # matrix = [[1, 0, -6, 0, -6, 3],
    #           [0, 1, 8, -1, 0, 2],
    #           [0, 0, 0, 0, 1, 0],
    #           [0, 0, 0, 0, 0, 0],
    #           [1, 0, -5, 0, -6, 3],
    #           [0, 1, 8, -1, 0, 2],
    #           [0, 0, 0, 0, 1, 0],
    #           [0, 0, 0, 0, 0, 0],
    #           [1, 0, -6, 0, -6, 3],
    #           [0, 1, 8, -1, 0, 2],
    #           [2, 0, 0, 0, 1, 0],
    #           [0, 0, 0, 0, 0, 0],
    #           [1, 0, -5, 0, -6, 3],
    #           [0, 1, 8, -1, 0, 2],
    #           [3, 0, 0, 0, 1, 0],
    #           [1, 0, -6, 0, -6, 3],
    #           [0, 1, 8, -1, 0, 2],
    #           [0, 0, 0, 0, 1, 0],
    #           [8, 0, 0, 0, 0, 0],
    #           [1, 0, -5, 0, -6, 3],
    #           [-3, 1, 8, -1, 0, 2],
    #           [0, 0, 0, 0, 1, 0]
    #           ]

    print("Starting matrix:")
    print_matrix(matrix)

    RREF: bool = True  # Reduced row echelon form (RREF) option

    init(autoreset=True) # Initialize colorama for error messages and auto-resetting after each line

    # Searches for the first nonzero entry in the first row to create a diagonal of expected pivots, checks for consistency and infinite number of solutions
    if not find_pivot_locations_and_store(matrix):
        print(Fore.RED + "FAILED: Failed to either locate pivots or verify the consistency of the augmented matrix.")
        return

    print(f"PIVOT COORDINATES: {pivot_coordinates}")
    expected_pivot_coordinates: List[Tuple[int, int, bool, bool]] = []

    row_index = 0
    column_index = 0
    #  To populate the expected pivot coordinates without considering the first nonzero entry for possible row exchange
    while row_index < (len(matrix)) and column_index < (len(matrix[0]) - 1):
            expected_location_value_equal_to_zero: bool = True if matrix[row_index][column_index] == 0 else False
            expected_location_value_equal_to_one: bool = True if matrix[row_index][column_index] == 1 else False
            expected_pivot_coordinates.append((row_index, column_index, expected_location_value_equal_to_one, expected_location_value_equal_to_zero))
            row_index, column_index = row_index + 1, column_index + 1
    print(f"EXPECTED PIVOT COORDINATES: {expected_pivot_coordinates}")

    # To check if the original pivot coordinates don't match the expected coordinates
    if not (pivot_coordinates == expected_pivot_coordinates):
        expected_coordinate_size = len(expected_pivot_coordinates)
        actual_coordinate_size = len(pivot_coordinates)
        if expected_coordinate_size - actual_coordinate_size > 0:
            # Ensure there exists an equal number of pivot coordinates as the expected to then be replaced if needed
            for i in range(actual_coordinate_size - 1, expected_coordinate_size - 1):
                pivot_coordinates.append(expected_pivot_coordinates[i + 1])

        for index, (expected_row, expected_column, expected_location_value_equal_to_one, expected_location_value_equal_to_zero) in enumerate(expected_pivot_coordinates):
            for i in range(expected_row + 1, len(matrix)):
                if expected_location_value_equal_to_zero and matrix[i][expected_column] != 0:
                    # Perform row exchange if a row below at the same column has any nonzero and the value at the expected pivot location is 0

                    # Exchange the expected row for the ith row if there's a 1 found
                    row_exchange(matrix, expected_row, i)

                    print(f"Matrix after row exchange:")
                    print_matrix(matrix)

                    # Find the first nonzero entry in the first row to create a new diagonal of pivots while verifying the consistency again
                    if not find_pivot_locations_and_store(matrix):
                        print(Fore.RED + "\nFAILED: Failed to either locate pivots or verify the consistency of the augmented matrix.")
                        return
                    break

    print(f"PIVOT COORDINATES AFTER EXCHANGE: {pivot_coordinates}\n")

    # To scale at expected pivot coordinates for nonones and nonzeros and eliminate nonzeros above and below pivots through row replacement
    while not pivots_equal_to_one(matrix):
        for (pivot_row, pivot_column, pivot_equal_to_one, pivot_equal_to_zero) in pivot_coordinates:
            # To create a pivot equal to 1 at each coordinate such that on the following row to the right, the element is also equal to 1 for each row that's not the last one
            pivot_value = matrix[pivot_row][pivot_column]
            pivot_equal_to_one, pivot_equal_to_zero = matrix[pivot_row][pivot_column] == 1, matrix[pivot_row][pivot_column] == 0

            if pivot_equal_to_zero:
                continue
            elif not pivot_equal_to_one:
                row_scaling(matrix, pivot_row, 1 / pivot_value, pivot_value)  # To change the pivot to 1

                print(f"Matrix after scaling:")
                print_matrix(matrix)

            # Search for nonzeros below each pivot and eliminate with pivot
            for i in range(pivot_row + 1, len(matrix)):
                # Find any nonzeros below the pivot to eliminate through row replacement
                value = matrix[i][pivot_column]
                if value != 0:
                    row_replacement(matrix, i, -value, pivot_row)
                else:
                   continue

            if pivot_row > 0 and RREF:
                for i in range(pivot_row - 1, -1, -1):
                    # Find any nonzeros above the pivot to eliminate through row replacement
                    value = matrix[i][pivot_column]
                    if value != 0:
                        row_replacement(matrix, i, -value, pivot_row)
                    else:
                        continue
                print()

            print("Matrix after row replacement:")
            print_matrix(matrix)

            if not find_pivot_locations_and_store(matrix):
                print(Fore.RED + "\nFAILED: Failed to either locate pivots or verify the consistency of the augmented matrix.")
                return

    if not infinite_number_of_solutions_validation(matrix):
        print(Fore.RED + "\nFAILED: An infinite number of solutions has been detected in the augmented matrix.")
        return

    if not RREF:
        print(Fore.GREEN + "SUCCESS: Transformed augmented matrix into row echelon form (REF).")
    else:
        print(Fore.GREEN + "SUCCESS: Transformed augmented matrix into reduced row echelon form (RREF).")

    print(Fore.GREEN + "\nFinal matrix:")
    print_matrix(matrix)


def row_scaling(matrix: List[List[float]], row_index: int, num: float, pivot_value: float = None) -> None: # The pivot is optional and only for printing the multiplier
    for index, value in enumerate(matrix[row_index]):
        matrix[row_index][index] = value * num

    if pivot_value is not None:
        print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by 1/{pivot_value if 0 < pivot_value < 1 else int(pivot_value)}.")
    else:
        print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by {num}.")

def row_replacement(matrix: List[List[float]], row_index: int, multiplier: float, pivot_row: int) -> None:
    if multiplier == 0:
        print(Fore.RED + f"ERROR: Cannot perform row replacement with a multiplier equal to 0!\n")

    for index, value in enumerate(matrix[row_index]):
        matrix[row_index][index] = value + multiplier * matrix[pivot_row][index] # Current value + multiplier * value above in pivot row

    print(Fore.GREEN + f"SUCCESS: Replaced row {row_index} with R{row_index} + {multiplier}*R{pivot_row}.")

def row_exchange(matrix: List[List[float]], row_index: int, second_row_index: int) -> None:
    temp_row: List[float] = matrix[row_index]
    matrix[row_index] = matrix[second_row_index]
    matrix[second_row_index] = temp_row

    print(Fore.GREEN + f"SUCCESS: Performed row exchange such that R{row_index} = R{second_row_index} and R{second_row_index} = R{row_index}_prev\n")

def find_pivot_locations_and_store(matrix: List[List[float]]) -> bool:
    global pivot_coordinates
    # To check if there's at least one nonzero entry in each row excluding the last element, if the matrix is consistent, and saves the first potential pivot
    non_zero_index_matrix: List[Tuple[int, int, bool, bool]] = []
    pivot_coordinates = []
    nonzero_rows = []

    # Populate the nonzero matrix with the locations of the nonzeros to find the first pivot and nonzero rows for infinite solutions check
    for row_index, row_value in enumerate(matrix):
        for col_index, col_value in enumerate(row_value):
            if col_value != 0 and row_index not in nonzero_rows:
                nonzero_equal_to_one, nonzero_equal_to_zero = col_value == 1, col_value == 0
                non_zero_index_matrix.append((row_index, col_index, nonzero_equal_to_one, nonzero_equal_to_zero))
                nonzero_rows.append(row_index)

    # Verify the consistency of the matrix
    for i in range(len(matrix)):
        if i not in nonzero_rows and matrix[i][-1] != 0:
            print(Fore.RED + "ERROR: No nonzero coefficient found, yet there's still a nonzero constant. This matrix is inconsistent.")
            return False
    print(Fore.GREEN + "SUCCESS: The augmented matrix is consistent. Proceeding with creating pivots and performing row operations if needed.")

    # Create a sorted nonzero coordinate list by column indices
    non_zero_index_matrix = sorted(non_zero_index_matrix, key=lambda tup: tup[1])
    first_pivot_coordinate = non_zero_index_matrix[0]

    # To populate the pivot coordinates list with the coordinates of each pivot
    if first_pivot_coordinate is not None:
        new_pivot_row, new_pivot_column, new_pivot_equal_to_one, new_pivot_equal_to_zero = first_pivot_coordinate
        prev_pivot_row, prev_pivot_column = new_pivot_row, new_pivot_column
        while (new_pivot_row != len(matrix)) and (new_pivot_column != len(matrix[new_pivot_row]) - 1):  # Start at the row of the first pivot and end at the last row of matrix
            if matrix[new_pivot_row][new_pivot_column] == 0 and new_pivot_column != (len(matrix[new_pivot_row]) - 2):
                new_pivot_column += 1
            elif matrix[new_pivot_row][new_pivot_column] == 0 and new_pivot_column == (len(matrix[new_pivot_row]) -2) and (prev_pivot_row + 1) != len(matrix):
                new_pivot_row += 1
                new_pivot_column = prev_pivot_column
            elif new_pivot_row != (len(matrix)) and new_pivot_column != (len(matrix[new_pivot_row]) - 1):
                new_pivot_equal_to_one = True if matrix[new_pivot_row][new_pivot_column] == 1 else False
                new_pivot_equal_to_zero = True if matrix[new_pivot_row][new_pivot_column] == 0 else False

                pivot_coordinates.append((new_pivot_row, new_pivot_column, new_pivot_equal_to_one, new_pivot_equal_to_zero))
                new_pivot_row, new_pivot_column = new_pivot_row + 1, new_pivot_column + 1
                prev_pivot_column = new_pivot_column
    else:
        print(Fore.YELLOW + "WARNING: No pivot found in the matrix.")

    print(Fore.GREEN + "\nSUCCESS: Located the pivots and verified the consistency of the augmented matrix.")
    return True

def infinite_number_of_solutions_validation(matrix: List[List[float]]) -> bool:
    # Search for columns only containing zeros, indicating an infinite number of solutions for the matrix
    for column in zip(*matrix):
        nonzeros_found = any(column_value != 0 for column_value in column)
        if nonzeros_found:  # Move onto the next column to search for columns only containing zeros
            continue
        print(Fore.RED + "ERROR: The augmented matrix has an infinite number of solutions.")
        return False

    print(Fore.GREEN + "SUCCESS: No infinite number of solutions detected in the augmented matrix.")
    return True

def pivots_equal_to_one(matrix: List[List[float]]) -> bool:
    # Search for nonzeros below each pivot and eliminate with pivot
    for (pivot_row, pivot_column, pivot_equal_to_one, pivot_equal_to_zero) in pivot_coordinates:
        if matrix[pivot_row][pivot_column] != 1:
            return False
    return True

def print_matrix(matrix: List[List[float]]):
    for row in matrix:
        print(row)

if __name__ == "__main__":
    main()