from quantum_circuit import QuantumCircuit as QC

if __name__ == "__main__":
    '''
    an example of grover's search with 3 qubits
    This link shows how it works: https://qiskit.org/textbook/ch-algorithms/grover.html
    '''
    grove_qc = QC(3)
    grove_qc.apply_hardmard(0)
    grove_qc.apply_hardmard(1)
    grove_qc.apply_hardmard(2)

    # oracle
    grove_qc.apply_controlZ(2, 0)
    grove_qc.apply_controlZ(1, 0)

    # amplification
    grove_qc.apply_hardmard(0)
    grove_qc.apply_hardmard(1)
    grove_qc.apply_hardmard(2)
    grove_qc.apply_pauliX(0)
    grove_qc.apply_pauliX(1)
    grove_qc.apply_pauliX(2)
    grove_qc.apply_controlZ([1, 2], 0)
    grove_qc.apply_pauliX(0)
    grove_qc.apply_pauliX(1)
    grove_qc.apply_pauliX(2)
    grove_qc.apply_hardmard(0)
    grove_qc.apply_hardmard(1)
    grove_qc.apply_hardmard(2)

    grove_qc.show_state()
    grove_qc.plot_pr()
