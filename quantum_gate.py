#!/usr/bin/env python
# coding: utf-8

import numpy as np

I = np.eye(2)

PX = np.array([[0, 1], 
               [1, 0]])

PY = np.array([[0, -np.complex(0, 1)], 
               [np.complex(0, 1), 0]])

PZ = np.array([[1, 0], 
               [0, -1]])

H = np.array([[1/np.sqrt(2), 1/np.sqrt(2)], 
              [1/np.sqrt(2), -1/np.sqrt(2)]])

SWAP = np.array([[1, 0, 0, 0], 
                [0, 0, 1, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1]])
