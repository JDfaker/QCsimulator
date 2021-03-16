from matplotlib import pyplot as plt
import numpy as np

from algorithm import Shors


def main(shors: Shors):
    '''

    Returns:

    '''
    thing = shors.compute_r()
    r = shors.cf(thing)
    while r % 2 != 0:
        r = shors.cf(shors.compute_r())
        a1 = int((shors.a ** (r / 2)) - 1)
        a2 = int((shors.a ** (r / 2)) + 1)
    a1 = int((shors.a ** (r / 2)) - 1)
    a2 = int((shors.a ** (r / 2)) + 1)
    return np.gcd(a1, shors.N), np.gcd(a2, shors.N)


if __name__ == "__main__":
    factors = []
    n = 78
    for i in range(0, 500):
        random = np.random.randint(0, n)
        while np.gcd(random, n) != 1:
            random = np.random.randint(0, n)

        shors = Shors(random, n)
        factors.append(main(shors))
    factors = np.array(factors).flatten()

    plt.hist(factors, bins=range(0, n + 1, 1))
    plt.xlabel('Factors')
    plt.ylabel('Frequencies')
    plt.title('Histogram of Shors Algorithm Output for  N = ' + str(n))
    plt.show()

