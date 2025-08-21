##########
#  Imports
import numpy as np

from qiskit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.synthesis import LieTrotter

##########################
# Constant Pauli Operators
X = SparsePauliOp("X")
Z = SparsePauliOp("Z")
Id = SparsePauliOp("I")

#######################
# Functions Definitions

def _calculate_H0(w00 : float, w01 : float, Axz : float) -> SparsePauliOp:
    """
    Calculates the time-independent Hamiltonian H0 for the provided parameters

    Parameters
    ----------
    w00 : float
        Larmor frequency of central qubit
    w01 : float
        Larmor frequency of coupled qubit
    Axz : float
        Hyperfine coupling between Sx and Iz

    Returns
    -------
    H0 : SparsePauliOp
        Time-independent Hamiltonian
    """
    H0 = w00/2*(Z ^ Id) + w01/2*(Id ^ Z) + Axz/4*(Z ^ X)
    return 2*np.pi*H0

def _calculate_H1(t : float, wp : float, w1 : float, phi : float)  -> SparsePauliOp:
    """
    Calculates the time-dependent Hamiltonian H1(t) during a pulse for the provided parameters

    Parameters
    ----------
    t : float
        Time variable
    wp : float
        Frequency of the pulse
    w1 : float
        Rabi frequency of the pulse
    phi : float
        Phase of the pulse

    Returns
    -------
    H1 : SparsePauliOp
        Time-dependent Hamiltonian
    """
    return 2*np.pi*w1*np.sin(2*np.pi*wp*t + phi)*(X ^ Id)

def _append_pulse(
        protocol : QuantumCircuit,
        duration : float,
        t0 : float,
        H0 : SparsePauliOp,
        wp : float,
        w1 : float,
        phi : float,
        dt : float
        ) -> None:
    """
    Appends a pulse to the protocol by applying several Pauli evolutions in sequence

    Parameters
    ----------
    protocol : QuantumCircuit
        Quantum circuit containing the protocol which the pulse will be appended
    duration : float
        Duration of the pulse
    t0 : float
        Initial time of the pulse
    H0 : SparsePauliOp
        Time-independent Hamiltonian
    wp : float
        Frequency of the pulse 
    w1 : float
        Rabi frequency of the pulse
    phi : float
        Phase of the pulse
    dt : float
        Time step for the time array of the pulse
    """
    tp_array = np.arange(t0, t0+duration, dt)
    for idx_tp in range(tp_array.size):
        protocol.append(PauliEvolutionGate(H0 + _calculate_H1(tp_array[idx_tp], wp, w1, phi),
                                           time=dt,
                                           synthesis=LieTrotter(reps=1)),
                                           [0,1])

def CPMG(
        tau : float,
        N : int,
        w00 : float,
        w01 : float,
        Axz : float,
        w1 : float,
        dt : float,
        wp=None,
        tpi=None
        ) -> QuantumCircuit:
    """
    Performs a CPMG sequences with the provided parameters by calling _append_pulse and Pauli evolution gates for the free evolutions.
    The CPMG is composed by pi_x pulses repeated N times, with an initial and final pi_y/2 pulse. 

    Parameters
    ----------
    N : int
        Number of pulses in the sequence
    w00 : float
        Larmor frequency of central qubit
    w01 : float
        Larmor frequency of coupled qubit
    Axz : float
        Hyperfine coupling between Sx and Iz
    w1 : float
        Rabi frequency of the pulse
    dt : float
        Time step for the time array of the pulse
    wp : None, float
        Frequency of the pulse. If None, assumes w00
    tpi : None, float
        Pulse duration. If None, assumes 1/(2*w1)

    Returns
    -------
    protocol : QuantumCircuit
        quantum circuit containing the protocol
    """
    if tpi is None:
        tpi = 1/(2*w1)
    if wp is None:
        wp = w00

    protocol = QuantumCircuit(2)
    H0 = _calculate_H0(w00, w01, Axz)

    # Pulse separation between pulses
    ps = tau - tpi
    if ps < 0:
        raise ValueError(f"Pulse separation cannot be negative. ps = tau - tpi = {tau - tpi}")

    # RY(pi/2)
    _append_pulse(protocol, tpi/2, 0, H0, wp, w1, -np.pi/2, dt)
    t0 = tpi/2

    # U(tau/2)
    protocol.append(PauliEvolutionGate(H0,
                                       time=ps/2-tpi/2,
                                       synthesis=LieTrotter(reps=1)),
                                       [0,1])
    t0 += ps/2-tpi/2

    for idx_N in range(N):
        # RX(pi)
        _append_pulse(protocol, tpi, t0, H0, wp, w1, 0, dt)
        t0 += tpi

        if idx_N != N-1:
            # U(tau)
            protocol.append(PauliEvolutionGate(H0,
                                               time=ps,
                                               synthesis=LieTrotter(reps=1)),
                                               [0,1])
            t0 += ps

    # U(tau/2)
    protocol.append(PauliEvolutionGate(H0,
                                       time=ps/2-tpi/2,
                                       synthesis=LieTrotter(reps=1)),
                                       [0,1])
    t0 += ps/2-tpi/2
        
    # RY(pi/2)
    _append_pulse(protocol, tpi/2, t0, H0, wp, w1, -np.pi/2, dt)

    return protocol

def XYN(
        tau : float,
        N : int,
        w00 : float,
        w01 : float,
        Axz : float,
        w1 : float,
        dt : float,
        wp = None,
        tpi=None
        ) -> QuantumCircuit:
    """
    Performs a XYN sequences with the provided parameters by calling _append_pulse and Pauli evolution gates for the free evolutions
    The XYN is composed by N intercalated pi_x and pi_y pulses, with an initial and final pi_y/2 pulse. 

    Parameters
    ----------
    N : int
        Number of pulses in the sequence
    w00 : float
        Larmor frequency of central qubit
    w01 : float
        Larmor frequency of coupled qubit
    Axz : float
        Hyperfine coupling between Sx and Iz
    w1 : float
        Rabi frequency of the pulse
    dt : float
        Time step for the time array of the pulse
    wp : None, float
        Frequency of the pulse. If None, assumes w00
    tpi : None, float
        Pulse duration. If None, assumes 1/(2*w1)

    Returns
    -------
    protocol : QuantumCircuit
        quantum circuit containing the protocol
    """
    if tpi is None:
        tpi = 1/(2*w1)
    if wp is None:
        wp = w00
    if N % 2 != 0:
        raise ValueError('N must be even in the XYN sequence')
    else:
        M = int(N/2)

    protocol = QuantumCircuit(2)
    H0 = _calculate_H0(w00, w01, Axz)

    # Pulse separation between pulses
    ps = tau - tpi
    if ps < 0:
        raise ValueError(f"Pulse separation cannot be negative. ps = tau - tpi = {tau - tpi}")

    # RY(pi/2)
    _append_pulse(protocol, tpi/2, 0, H0, wp, w1, -np.pi/2, dt)
    t0 = tpi/2

    # U(tau/2)
    protocol.append(PauliEvolutionGate(H0,
                                       time=ps/2-tpi/2,
                                       synthesis=LieTrotter(reps=1)),
                                       [0,1])
    t0 += ps/2-tpi/2

    for idx_M in range(M):
        # RX(pi)
        _append_pulse(protocol, tpi, t0, H0, wp, w1, 0, dt)
        # U(tau)
        protocol.append(PauliEvolutionGate(H0,
                                            time=ps,
                                            synthesis=LieTrotter(reps=1)),
                                            [0,1])
        t0 += ps

        # RY(pi)
        _append_pulse(protocol, tpi, t0, H0, wp, w1, -np.pi/2, dt)
        t0 += tpi

        if idx_M != N-1:
            # U(tau)
            protocol.append(PauliEvolutionGate(H0,
                                               time=ps,
                                               synthesis=LieTrotter(reps=1)),
                                               [0,1])
            t0 += ps

    # U(tau/2)
    protocol.append(PauliEvolutionGate(H0,
                                       time=ps/2-tpi/2,
                                       synthesis=LieTrotter(reps=1)),
                                       [0,1])
    t0 += ps/2-tpi/2
        
    # RY(pi/2)
    _append_pulse(protocol, tpi/2, t0, H0, wp, w1, -np.pi/2, dt)

    return protocol