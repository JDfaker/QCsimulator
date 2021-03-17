#!/usr/bin/env python
# coding: utf-8

from quantum_circuit import QuantumCircuit as QC

'''
an example of grove search with 3 qubit
This link show how it works: https://qiskit.org/textbook/ch-algorithms/grover.html
'''
if __name__ == "__main__":
    my_qc = QC(5)
    my_qc.apply_hardmard(0)
    my_qc.apply_hardmard(1)
    my_qc.apply_hardmard(2)
    my_qc.apply_hardmard(3)
    my_qc.apply_hardmard(4)

    my_qc.apply_grover_oracle(2)
    my_qc.apply_amplification()
    my_qc.apply_grover_oracle(2)
    my_qc.apply_amplification()
    my_qc.apply_grover_oracle(2)
    my_qc.apply_amplification()
    my_qc.apply_grover_oracle(2)
    my_qc.apply_amplification()
    my_qc.apply_grover_oracle(2)
    my_qc.apply_amplification()
    my_qc.plot_pr()
