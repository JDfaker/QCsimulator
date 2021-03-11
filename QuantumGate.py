#!/usr/bin/env python
# coding: utf-8

import numpy as np
from sparse_matrix import SparseMatrix

I = SparseMatrix.sparsify(np.eye(2))

PX = SparseMatrix.sparsify(np.array([[0, 1],
               [1, 0]]))

PY = SparseMatrix.sparsify((np.array([[0, -np.complex(0, 1)],
               [np.complex(0, 1), 0]]))

PZ = SparseMatrix.sparsify(np.array([[1, 0],
               [0, -1]]))

H = SparseMatrix.sparsify(np.array([[1 / np.sqrt(2), 1 / np.sqrt(2)],
              [1 / np.sqrt(2), -1 / np.sqrt(2)]]))

SWAP = SparseMatrix.sparsify(np.array([[1, 0, 0, 0],
                 [0, 0, 1, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1]]))
