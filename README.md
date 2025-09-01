## **Gaussian Elimination Implementation by Pedro Palacios for Linear Algebra (MAT-2233) at UT San Antonio**

This program was completed as extra-credit to demonstrate an understanding of elementary row operations on an augmented matrix to transform it into either row echelon form (REF)
or reduced row echelon form (RREF). The augmented matrix is set in the program, however, it can be configured to retrieve user input; previously, I used numpy's rng integers method
to generate a 4 x 4 matrix with random integers. 

Currently, this program has the following row operations:
  1. Row Scaling
  2. Row Replacement

Included are error messages and warnings that will output if the matrix is inconsistent (coefficients of zero on variables equaling a nonzero), no pivots are found, or if
a value is being multiplied by a value of zero during row replacement.

