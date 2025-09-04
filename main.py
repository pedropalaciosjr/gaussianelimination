# Gaussian elimination implementation by Pedro Palacios for MAT-2233 at UT San Antonio
from colorama import Fore, init
from typing import Tuple, List

def main() -> None:
    global pivot_coordinates
    matrix: List[List[float]] = [[0, 1, 5, 3],
                                [1, 4, 4, -3],
                                [2, 6, 3, -2]]
    print("Starting matrix:")
    print_matrix(matrix)

    RREF: bool = True  # Reduced row echelon form (RREF) option

    init(autoreset=True) # Initialize colorama for error messages and auto-resetting after each line

    # Searches for the first nonzero entry in the first row to create a diagonal of expected pivots and checks for consistency
    find_pivot_locations_and_store(matrix)

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
                    find_pivot_locations_and_store(matrix)
                    break

    print(f"PIVOT COORDINATES AFTER EXCHANGE: {pivot_coordinates}\n")

    # To scale at expected pivot coordinates for nonones and nonzeros and eliminate nonzeros above and below pivots through row replacement
    for (pivot_row, pivot_column, pivot_equal_to_one, pivot_equal_to_zero) in pivot_coordinates:
        # To create a pivot equal to 1 at each coordinate such that on the following row to the right, the element is also equal to 1 for each row that's not the last one
        pivot_value = matrix[pivot_row][pivot_column]

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
        print(Fore.RED + "\nFAILED: Failed to locate pivots or verify the consistency of the augmented matrix.")
        return

    print(Fore.GREEN + "\nSUCCESS: Located the pivots and verified the consistency of the augmented matrix.")
    print(Fore.GREEN + "SUCCESS: Transformed augmented matrix into reduced row echelon form (RREF).")
    print(Fore.YELLOW + "\nFinal matrix:")
    print_matrix(matrix)

def row_scaling(matrix: List[List[float]], row_index: int, num: float, pivot_value: float = None) -> None: # The pivot is optional and only for printing the multiplier
    for index, value in enumerate(matrix[row_index]):
        matrix[row_index][index] = value * num

    if pivot_value is not None:
        print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by 1/{pivot_value if 0 < pivot_value < 1 else int(pivot_value)}.")
    else:
        print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by {num}")

def row_replacement(matrix: List[List[float]], row_index: int, multiplier: float, pivot_row: int) -> None:
    if multiplier == 0:
        print(Fore.RED + f"ERROR: Cannot perform row replacement with a multiplier equal to 0!\n")

    for index, value in enumerate(matrix[row_index]):
        matrix[row_index][index] = value + multiplier * matrix[pivot_row][index] # Current value + multiplier * value above in pivot row

    print(Fore.GREEN + f"SUCCESS: Replaced row {row_index} with R{row_index} + {multiplier}*R{pivot_row}")

def row_exchange(matrix: List[List[float]], row_index: int, second_row_index: int) -> None:
    temp_row: List[float] = matrix[row_index]
    matrix[row_index] = matrix[second_row_index]
    matrix[second_row_index] = temp_row

    print(Fore.GREEN + f"SUCCESS: Performed row exchange such that R{row_index} = R{second_row_index} and R{second_row_index} = R{row_index}_prev\n")

def find_pivot_locations_and_store(matrix: List[List[float]]) -> bool:
    global pivot_coordinates
    # To check if there's at least one nonzero entry in each row excluding the last element, if the matrix is consistent, and saves the first potential pivot
    non_zero_found = False
    first_pivot_coordinate = None
    nonzero_indicies: List[Tuple[int, int, bool, bool]] = []
    non_zero_index_matrix: List[Tuple[int, int]] = [(row_index, col_index) if col_value != 0 else (-1, -1) for row_index, row_value in enumerate(matrix) for col_index, col_value in enumerate(row_value)]
    pivot_coordinates = []

    for i in range(len(matrix)):
        for (row, col) in non_zero_index_matrix:
            if row == i and col != len(matrix[i]):  # If there exists a nonzero coordinate in the ith row that's not the last element
                non_zero_found = True
                equal_to_one = matrix[row][col] == 1
                equal_to_zero = matrix[row][col] == 0

                # To add the coordinates of potential pivots, the first nonzero entry of each row
                nonzero_indicies.append((row, col, equal_to_one, equal_to_zero))

                first_pivot_coordinate = nonzero_indicies[0]
        if not non_zero_found and matrix[i][-1] != 0:
            # If there's only zeros in the coefficients and a nonzero in the row's last column—inconsistent
            print(Fore.RED + "ERROR: No nonzero coefficient found, yet there's still a nonzero constant. This matrix is inconsistent.")
            return False
    print(Fore.GREEN + "SUCCESS: The augmented matrix is consistent. Proceeding with creating pivots and performing row operations for RREF.")

    # To populate the pivot coordinates list with the coordinates of each pivot
    if first_pivot_coordinate is not None:
        new_pivot_row, new_pivot_column, new_pivot_equal_to_one, new_pivot_equal_to_zero = first_pivot_coordinate
        for i in range(first_pivot_coordinate[0], len(matrix)):  # Start at the row of the first pivot and end at the last row of matrix
            if matrix[new_pivot_row][new_pivot_column] == 0 and new_pivot_column != (len(matrix[new_pivot_row]) - 2):
                new_pivot_column += 1
            elif new_pivot_row != (len(matrix)) and new_pivot_column != (len(matrix[new_pivot_row]) - 1):
                new_pivot_equal_to_one = True if matrix[new_pivot_row][new_pivot_column] == 1 else False
                new_pivot_equal_to_zero = True if matrix[new_pivot_row][new_pivot_column] == 0 else False
                pivot_coordinates.append((new_pivot_row, new_pivot_column, new_pivot_equal_to_one, new_pivot_equal_to_zero))
                new_pivot_row, new_pivot_column = new_pivot_row + 1, new_pivot_column + 1
    else:
        print(Fore.YELLOW + "WARNING: No pivot found in the matrix.")

    return True

def print_matrix(matrix: List[List[float]]):
    for row in matrix:
        print(row)

if __name__ == "__main__":
    main()