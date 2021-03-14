"""
This module holds the Sparse Matrix class
"""

from __future__ import annotations
from typing import List
import numpy as np


class SparseMatrix:
    """
    This is a class for sparse matrix interface. \n
    It holds the values in an inner numpy array in the following format: \n

    r_1 | r_2 | r_3 | ... | rn \n
    c_1 | c_2 | c_3 | ... | cn \n
    v_1 | v_2 | v_3 | ... | vn \n

    where r is the row index, c - column index and v is the values stored.

    The class implements the matrix mathematics for sparse matrices, \n
    also it provides conversion tools for conversion
    between numpy array and Sparse Matrix.
    """

    def __init__(self, values: List, rows: int, cols: int):
        """
        This is the constructor method for the Sparse Matrix \n
        @param values: List of tuples of (row, column, val) to specify the row
        and column where the val v is inserted \n
        @param rows: Number of rows in the matrix \n
        @param cols: Number of rows in the matrix \n
        """
        self.rows = rows
        self.cols = cols
        self.inner_array = np.zeros((3, len(values)))
        indexes = set()
        # Sort by row and then by column
        for i, (row, col, val) in enumerate(
                sorted(values, key=lambda x: (x[0], x[1]))):
            # Check that there are no duplicate entries for rth row and cth col
            assert (row, col) not in indexes, \
                "Found duplicate entry for row {}, column {}".format(row, col)

            indexes.add((row, col))
            self.inner_array[0][i] = row
            self.inner_array[1][i] = col
            self.inner_array[2][i] = val

    @staticmethod
    def sparsify(matrix: np.array) -> SparseMatrix:
        """
        Static method for generating Sparse Matrix for numpy matrix \n
        @param matrix: Numpy matrix \n
        @return: Sparse matrix generated \n
        """
        values = []
        rows, cols = np.nonzero(matrix)
        for i, _ in enumerate(rows):
            row, col = int(rows[i]), int(cols[i])
            values.append((row, col, matrix[row][col]))

        rows, cols = matrix.shape
        return SparseMatrix(values, rows, cols)

    def numpy(self) -> np.array:
        """
        Convert Sparse Matrix to numpy array \n
        @return: Numpy array \n
        """
        numpy_matrix = np.zeros((self.rows, self.cols))
        for i in range(len(self.inner_array[0])):
            row = int(self.inner_array[0][i])
            col = int(self.inner_array[1][i])
            val = self.inner_array[2][i]
            numpy_matrix[row][col] = val
        return numpy_matrix

    def tensordot(self, matrix: SparseMatrix) -> SparseMatrix:
        """
        Tensor product method for Sparse Matrices \n
        @param matrix: Sparse Matrix to be dotted with \n
        @return: Tensor product of the Sparse Matrices \n
        """
        values = []
        for i in range(self.inner_array.shape[1]):
            row_a = int(self.inner_array[0][i])
            col_a = int(self.inner_array[1][i])
            val_a = self.inner_array[2][i]
            for j in range(matrix.inner_array.shape[1]):
                row_b = int(matrix.inner_array[0][j])
                col_b = int(matrix.inner_array[1][j])
                val_b = matrix.inner_array[2][j]

                row_final = (row_a)*(matrix.rows) + row_b
                col_final = (col_a)*(matrix.cols) + col_b
                value = val_a * val_b

                values.append((int(row_final), int(col_final), value))

        rows, cols = self.rows*matrix.rows, self.cols*matrix.cols

        return SparseMatrix(values, rows, cols)

    def dot(self, matrix: SparseMatrix) -> SparseMatrix:
        """
        Dot product between two Sparse Matrices \n
        @param matrix: Sparse Matrix to be dotted with \n
        @return: Result Sparse Matrix \n
        """
        assert self.cols == matrix.rows, \
            "Dimensions do not match, {} != {}".format(self.cols, self.rows)

        values = []
        for row in self.get_nonzero_rows():
            row_vals = self.get_row(row)

            for col in matrix.get_nonzero_cols():
                col_vals = matrix.get_col(col)
                val = 0

                for c_val in row_vals:
                    if c_val in col_vals.keys():
                        val += row_vals[c_val] * col_vals[c_val]
                        if val != 0:
                            values.append((row, col, val))

        return SparseMatrix(values, self.rows, matrix.cols)

    def multiply(self, multiply: int) -> SparseMatrix:
        """
        Matrix multiplied by an integer for Sparse Matrices \n
        @param multiply: Integer to be multiplied with matrix \n
        @return: Product of the Sparse Matrix and integer \n
        """

        values = []
        for i in range(len(self.inner_array[2])):
            values.append((int(self.inner_array[0][i]),
                           int(self.inner_array[1][i]),
                           self.inner_array[2][i] * multiply))
        return SparseMatrix(values, self.rows, self.cols)

    def minus(self, matrix: SparseMatrix) -> SparseMatrix:
        """
        Method for subtracting Sparse Matrices \n
        @param matrix: Sparse Matrix to be subtracted \n
        @return: Result of subtracted matrices \n
        """
        a = SparseMatrix.numpy(self)
        b = SparseMatrix.numpy(matrix)

        return SparseMatrix.sparsify((a-b))

    def transpose(self) -> SparseMatrix:
        """
        Method for transposing a Sparse Matrices \n
        @return: Transposed sparse matrix \n
        """
        rows = self.cols
        cols = self.rows

        values = []
        for i in range(len(self.inner_array[2])):
            row = self.inner_array[1][i]
            col = self.inner_array[0][i]
            val = self.inner_array[2][i]
            values.append((int(row), int(col), val))

        return SparseMatrix(values, rows, cols)

    def get_row(self, row: int) -> dict:
        """
        Get all values in the row r \n
        @param r: Selected row \n
        @return: Dictionary of values in format vals[col]
                    with value at row row and column col \n
        """
        vals = {}
        for i in range(self.inner_array.shape[1]):
            # The inner array is sorted by row and column
            # This reduces the number of iterations
            if row < self.inner_array[0][i]:
                break

            if self.inner_array[0][i] == row:
                vals[self.inner_array[1][i]] = self.inner_array[2][i]
        return vals

    def get_col(self, col: int) -> dict:
        """
        Get all values in the column c \n
        @param c: Selected column \n
        @return: Dictionary of values in format vals[r]
                    with value at row rrow and column col \n
        """
        vals = {}
        for i in range(self.inner_array.shape[1]):
            if self.inner_array[1][i] == col:
                vals[self.inner_array[0][i]] = self.inner_array[2][i]
        return vals

    def get_value(self, row: int, col: int) -> float:
        """
        @param r: Selected row \n
        @param c: Selected column \n
        @return: Value at selected row and column \n
        """
        for i in range(self.inner_array.shape[1]):
            # The inner array is sorted by row and column
            # This reduces the number of iterations
            if row < self.inner_array[0][i]:
                break

            elif (row == self.inner_array[0][i] and
                  col == self.inner_array[1][i]):
                return self.inner_array[2][i]
        return 0

    def get_nonzero_rows(self) -> List[float]:
        """
        Get all rows that have entries \n
        @return: List of rows \n
        """
        return np.unique(self.inner_array[0]).tolist()

    def get_nonzero_cols(self) -> List[float]:
        """
        Get all columns that have entries \n
        @return: List of all columns \n
        """
        return np.unique(self.inner_array[1]).tolist()
