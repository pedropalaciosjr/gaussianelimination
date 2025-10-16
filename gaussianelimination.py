# Gaussian elimination implementation by Pedro Palacios for MAT-2233 at UT San Antonio
from colorama import Fore, init
from typing import Tuple, List
import sys

class GaussianElimination:
    def __init__(self, matrix) -> None:
        self.matrix = matrix
        self.gaussianelimination()

    def gaussianelimination(self) -> None:
        global pivot_coordinates


        print("Starting matrix:")
        self.print_matrix()

        RREF: bool = True  # Reduced row echelon form (RREF) option

        init(autoreset=True) # Initialize colorama for error messages and auto-resetting after each line

        # Searches for nonzero entries in each row and sorts them by column such that the element closest to the left and highest will be the pivot
        if not self.find_pivot_locations_and_store():
            print(Fore.RED + "\nFAILED: Failed to either locate pivots or verify the consistency of the augmented matrix.")
            return

        print(f"PIVOT COORDINATES: {pivot_coordinates}")
        expected_pivot_coordinates: List[Tuple[int, int, bool, bool]] = []

        row_index = 0
        column_index = 0
        #  To populate the expected pivot coordinates without considering the first nonzero entry for possible row exchange
        while row_index < (len(self.matrix)) and column_index < (len(self.matrix[0]) - 1):
                expected_location_value_equal_to_zero: bool = True if self.matrix[row_index][column_index] == 0 else False
                expected_location_value_equal_to_one: bool = True if self.matrix[row_index][column_index] == 1 else False
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
                for i in range(expected_row + 1, len(self.matrix)):
                    if expected_location_value_equal_to_zero and self.matrix[i][expected_column] != 0:
                        # Perform row exchange if a row below at the same column has any nonzero and the value at the expected pivot location is 0

                        # Exchange the expected row for the ith row if there's a 1 found
                        self.row_exchange(expected_row, i)

                        print(f"Matrix after row exchange:")
                        self.print_matrix()

                        # Find the first nonzero entry in the first row to create a new diagonal of pivots while verifying the consistency again
                        if not self.find_pivot_locations_and_store():
                            print(Fore.RED + "\nFAILED: Failed to either locate pivots or verify the consistency of the augmented matrix.")
                            return
                        break

        print(f"PIVOT COORDINATES AFTER EXCHANGE: {pivot_coordinates}\n")

        # To scale at expected pivot coordinates for nonones and nonzeros and eliminate nonzeros above and below pivots through row replacement
        while not self.echelon_form(RREF):
            for (pivot_row, pivot_column, pivot_equal_to_one, pivot_equal_to_zero) in pivot_coordinates:
                # To create a pivot equal to 1 at each coordinate such that on the following row to the right, the element is also equal to 1 for each row that's not the last one
                pivot_value = self.matrix[pivot_row][pivot_column]
                pivot_equal_to_one, pivot_equal_to_zero = self.matrix[pivot_row][pivot_column] == 1, self.matrix[pivot_row][pivot_column] == 0

                if pivot_equal_to_zero:
                    continue
                elif not pivot_equal_to_one:
                    self.row_scaling(pivot_row, 1 / pivot_value, pivot_value)  # To change the pivot to 1

                    print(f"Matrix after scaling:")
                    self.print_matrix()

                # Search for nonzeros below each pivot and eliminate with pivot
                for i in range(pivot_row + 1, len(self.matrix)):
                    # Find any nonzeros below the pivot to eliminate through row replacement
                    value = self.matrix[i][pivot_column]
                    if value != 0:
                        self.row_replacement(i, -value, pivot_row)
                    else:
                       continue

                if pivot_row > 0 and RREF:
                    for i in range(pivot_row - 1, -1, -1):
                        # Find any nonzeros above the pivot to eliminate through row replacement
                        value = self.matrix[i][pivot_column]
                        if value != 0:
                            self.row_replacement(i, -value, pivot_row)
                        else:
                            continue
                    print()

                print("Matrix after row replacement:")
                self.print_matrix()

                if not self.find_pivot_locations_and_store():
                    print(Fore.RED + "\nFAILED: Failed to either locate pivots or verify the consistency of the augmented matrix.")
                    return

        if not self.infinite_number_of_solutions_validation():
            print(Fore.RED + "\nFAILED: An infinite number of solutions has been detected in the augmented matrix.")
            return

        if not RREF:
            print(Fore.GREEN + "SUCCESS: Transformed augmented matrix into row echelon form (REF).")
        else:
            print(Fore.GREEN + "SUCCESS: Transformed augmented matrix into reduced row echelon form (RREF).")

        print(Fore.GREEN + "\nFinal matrix:")
        self.print_matrix()


    def row_scaling(self, row_index: int, num: float, pivot_value: float = None) -> None: # The pivot is optional and only for printing the multiplier
        if num == 0:
            print(Fore.RED + f"ERROR: Cannot perform row scaling with a value equal to 0!\n")
            sys.exit()

        for index, value in enumerate(self.matrix[row_index]):
            self.matrix[row_index][index] = value * num

        if pivot_value is not None:
            print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by 1/{pivot_value if 0 < pivot_value < 1 else int(pivot_value)}.")
        else:
            print(Fore.GREEN + f"SUCCESS: Scaled row {row_index} by {num}.")

    def row_replacement(self, row_index: int, multiplier: float, pivot_row: int) -> None:
        if multiplier == 0:
            print(Fore.RED + f"ERROR: Cannot perform row replacement with a multiplier equal to 0!\n")
            sys.exit()

        for index, value in enumerate(self.matrix[row_index]):
            self.matrix[row_index][index] = value + multiplier * self.matrix[pivot_row][index] # Current value + multiplier * value above in pivot row

        print(Fore.GREEN + f"SUCCESS: Replaced row {row_index} with R{row_index} + {multiplier}*R{pivot_row}.")

    def row_exchange(self, row_index: int, second_row_index: int) -> None:
        temp_row: List[float] = self.matrix[row_index]
        self.matrix[row_index] = self.matrix[second_row_index]
        self.matrix[second_row_index] = temp_row

        print(Fore.GREEN + f"SUCCESS: Performed row exchange such that R{row_index} = R{second_row_index} and R{second_row_index} = R{row_index}_prev\n")

    def find_pivot_locations_and_store(self) -> bool:
        global pivot_coordinates
        # To check if there's at least one nonzero entry in each row excluding the last element, if the matrix is consistent, and saves the first potential pivot
        non_zero_index_matrix: List[Tuple[int, int, bool, bool]] = []
        pivot_coordinates = []
        nonzero_rows = []

        # Populate the nonzero matrix with the locations of the nonzeros to find the first pivot and nonzero rows for infinite solutions check
        for row_index, row_value in enumerate(self.matrix):
            for col_index, col_value in enumerate(row_value):
                if col_value != 0 and row_index not in nonzero_rows:
                    nonzero_equal_to_one, nonzero_equal_to_zero = col_value == 1, col_value == 0
                    non_zero_index_matrix.append((row_index, col_index, nonzero_equal_to_one, nonzero_equal_to_zero))
                    nonzero_rows.append(row_index)

        if len(non_zero_index_matrix) == 0:
            print(Fore.YELLOW + "WARNING: No nonzeros found in the matrix, so pivots couldn't be located.")
            return False

        # Verify the consistency of the matrix
        for i in range(len(self.matrix)):
            if i not in nonzero_rows and self.matrix[i][-1] != 0:
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
            while (new_pivot_row != len(self.matrix)) and (new_pivot_column != len(self.matrix[new_pivot_row]) - 1):  # Start at the row of the first pivot and end at the last row of matrix
                if self.matrix[new_pivot_row][new_pivot_column] == 0 and new_pivot_column != (len(self.matrix[new_pivot_row]) - 2):
                    new_pivot_column += 1
                elif self.matrix[new_pivot_row][new_pivot_column] == 0 and new_pivot_column == (len(self.matrix[new_pivot_row]) -2) and (prev_pivot_row + 1) != len(self.matrix):
                    new_pivot_row += 1
                    new_pivot_column = prev_pivot_column
                elif new_pivot_row != (len(self.matrix)) and new_pivot_column != (len(self.matrix[new_pivot_row]) - 1):
                    new_pivot_equal_to_one = True if self.matrix[new_pivot_row][new_pivot_column] == 1 else False
                    new_pivot_equal_to_zero = True if self.matrix[new_pivot_row][new_pivot_column] == 0 else False

                    pivot_coordinates.append((new_pivot_row, new_pivot_column, new_pivot_equal_to_one, new_pivot_equal_to_zero))
                    new_pivot_row, new_pivot_column = new_pivot_row + 1, new_pivot_column + 1
                    prev_pivot_column = new_pivot_column
        else:
            print(Fore.YELLOW + "WARNING: No pivot found in the matrix.")

        print(Fore.GREEN + "\nSUCCESS: Located the pivots and verified the consistency of the augmented matrix.")
        return True

    def infinite_number_of_solutions_validation(self) -> bool:
        # Search for columns only containing zeros, indicating an infinite number of solutions for the matrix
        for column in zip(*self.matrix):
            nonzeros_found = any(column_value != 0 for column_value in column)
            if nonzeros_found:  # Move onto the next column to search for columns only containing zeros
                continue
            print(Fore.RED + "ERROR: The augmented matrix has an infinite number of solutions.")
            return False

        print(Fore.GREEN + "SUCCESS: No infinite number of solutions detected in the augmented matrix.")
        return True

    def echelon_form(self, RREF: bool) -> bool:
        # Check if the matrix is in RREF or REF by checking the value of the pivots and values above and/or below
        for (pivot_row, pivot_column, pivot_equal_to_one, pivot_equal_to_zero) in pivot_coordinates:
            if self.matrix[pivot_row][pivot_column] != 1:
                return False
            for index, column in enumerate(zip(*self.matrix)):
                if index != pivot_column:
                    continue
                if not RREF:
                    column = column[pivot_row:]
                nonzeros = []
                for column_value in column:
                    if column_value != 0:
                        nonzeros.append(column_value)

                if len(nonzeros) > 1:
                    return False

        return True

    def print_matrix(self):
        for row in self.matrix:
            print(row)

    def get_matrix(self) -> List[List[float]]:
        return self.matrix