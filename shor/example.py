"""
A script for example use of classes Shors, ShorsSecondRegioster, and QFT \n
Outputs are graphs showing features of classes \n
"""

from matplotlib import pyplot as plt
import numpy as np
import shors_testing as st
import time
from algorithm import Shors

factors = []
n = 78
for i in range(0, 500):
    random = np.random.randint(0, n)
    while np.gcd(random, n) != 1:
        random = np.random.randint(0, n)

    shors = Shors(random, n)
    factors.append(Shors.compute_shors(shors))
factors = np.array(factors).flatten()

plt.hist(factors, bins=range(0, n + 1, 1))
plt.xlabel('Factors')
plt.ylabel('Frequencies')
plt.title('Histogram of Shors Algorithm Output for  N = ' + str(n))
plt.show()

st.run_time(8)
st.test_outputs(10)
