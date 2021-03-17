import numpy as np
import time
import sympy
from algorithm import Shors
from quantum_fourier_transform import QFT
from second_register import ShorsSecondRegister
import matplotlib.pyplot as plt


def run_time(qubit_no):
    """
    Function to calculate the run time of Shors \n
    @param quibit_no: Integer number of qubits \n
    @return: Graph of time for iterations \n
    """
    times = []
    for i in range(3, qubit_no):
        value = 2**i
        random = 0
        datapoints = 0
        start = time.process_time()
        random = np.random.randint(0, value)
        while np.gcd(random, i) != 1:
            random = np.random.randint(0, value)
        datapoints = np.array(Shors(random, value).compute_shors())
        times.append(time.process_time() - start)
    plt.figure(figsize=(10, 5))
    plt.plot(np.arange(3, (len(times)+3)), np.array(times))
    plt.ylabel('Time(s)')
    plt.xlabel('Qubit Number')
    plt.title('Time for a single itteration as a function of # qubits')
    plt.show()


def test_outputs(n):
    """
    Function to calculate number correct factors determined \n
    @param n: Number of qubits in the system \n
    @return: Graph of number of correct factors \n
    """
    success = []
    factors_true = []
    accuracy = []
    inaccuracy = []
    for i in range(0, n):
        value = np.fromiter(sympy.ntheory.factorint(i).keys(), dtype=int)
        factors_true.append(value.tolist())
    for i in range(5, n):
        factors_new = []
        success_temp = 0
        for j in range(0,50):
            random = np.random.randint(0, i)
            while np.gcd(random, n) != 1:
                random = np.random.randint(0, i)
            datapoints = np.array(Shors(random, i).compute_shors())
            factors_new.append(datapoints)
        factors_new = np.unique(np.array(factors_new).flatten())
        for k in range(len(factors_true[i])):
            if np.sum(np.where(factors_new % factors_true[i][k] == 0)) >= 1:
                success_temp += 1
            else:
                pass
        accuracy.append(success_temp/len(factors_new)*100)
        inaccuracy.append((1-(success_temp/len(factors_new)))*100)
    fig, ax = plt.subplots(figsize=(15, 5))
    labels = np.arange(5, n)
    ax.bar(labels, accuracy, width=0.5, label='Correct')
    ax.bar(labels, inaccuracy, width=0.5, bottom=accuracy,
           label='Incorrect')
    ax.legend()
    ax.set_title('Number of Correct Factors determinedfor $ 5 < N \leq 50$')
    ax.set_ylabel('Percentage of factors')
    ax.set_xlabel('N')
    plt.show()
