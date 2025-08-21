##########
#  Imports
import numpy as np
from dd_gates import CPMG

from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

##############
# IBMQ account
from qiskit_ibm_runtime import QiskitRuntimeService
service = QiskitRuntimeService(channel="ibm_quantum_platform",
                               # your token
                               token="", 
                               # your instance
                               instance='') 

#####################
# Transpilation Setup

real_backend = service.backend('ibm_torino')
estimator = Estimator(real_backend)
estimator.options.default_shots = 100
pass_manager = generate_preset_pass_manager(
     # Optimization with inverse cancelation but no commutative cancellation
    optimization_level=2,
    # Use the backend to transpile
    backend=real_backend, 
    # Take the qubit with longer T2 and smallest gate time
    layout_method='dense', 
    # Use the synthesis method to map the gates to the native operations of the backend,
    # translator method can also be used but results in longer pulses
    translation_method="synthesis", 
    # No approximation of the gates
    approximation_degree=1, 
                                        )

###########
# Variables

w00 = 50
w01 = 1
Axz = 0.1
w1 = 5
dt = 0.001

Z = SparsePauliOp("Z")
Id = SparsePauliOp("I")
observables = [Z^Id, Id^Z]

tau = np.arange(0.1, 4, 0.01)
N = 10

####################
# Generation of PUBs
pubs = []

# Generation of the PUBs, each one for a value of tau
for val_tau in tau:
    transpiled_protocol = pass_manager.run(
        CPMG(val_tau, N, w00, w01, Axz, w1, dt)
        )
    mapped_observables = [
        obs.apply_layout(transpiled_protocol.layout) for obs in observables
        ]
    pubs.append((transpiled_protocol, mapped_observables))

###############
# Job Execution
job = estimator.run(pubs)

print('Sent Job ID:', job.job_id())