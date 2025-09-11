## **Gaussian Elimination Implementation by Pedro Palacios for Linear Algebra (MAT-2233) at UT San Antonio**

This program was completed as extra credit to demonstrate an understanding of elementary row operations on an augmented matrix to transform it into either row echelon form (REF) or reduced row echelon form (RREF). The variable RREF can be manually set in the program to transform the matrix into RREF or REF. The augmented matrix is set in the program, however, it can be configured to retrieve user input, or NumPy can be used to randomly generate matrices. There are six provided augmented matrices; you can uncomment whichever one to test and comment the others.

Currently, this program has the following row operations:
  1. Row Scaling
  2. Row Replacement
  3. Row Exchange

Included are error messages and warnings that will output if the matrix is inconsistent (coefficients of zero on variables equaling a nonzero), no pivots are found, a value is being multiplied by a value of zero during row replacement, or an augmented matrix has an infinite number of solutions (defined such that the number of variables > number of nonzero rows). 

