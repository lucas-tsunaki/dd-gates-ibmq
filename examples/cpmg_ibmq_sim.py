##########
#  Imports
import numpy as np
import matplotlib.pyplot as plt
from dd_gates import CPMG

from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator

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
aer_backend = AerSimulator.from_backend(real_backend)
estimator = Estimator(aer_backend)
pass_manager = generate_preset_pass_manager(
     # Optimization with inverse cancelation but no commutative cancellation
    optimization_level=2,
    # Use the backend to transpile
    backend=aer_backend, 
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
print('Running')
job = estimator.run(pubs)
job_result = job.result()

######################
# Retrieve the Results
expectation_values = np.empty((len(tau), len(observables)), dtype=float)
stds = np.empty((len(tau), len(observables)), dtype=float)

for idx_tau, pub_result in enumerate(job_result):
    for idx_obs, _ in enumerate(observables):
        expectation_values[idx_tau, idx_obs] = pub_result.data.evs[idx_obs]
        stds[idx_tau, idx_obs] = pub_result.data.stds[idx_obs]

np.savetxt('./cpmg_ibmq_sim_results.txt', np.c_[tau, expectation_values, stds], header='tau\tobs1\tobs2\tstd1\tstd2')

##############
# Plot Results
plt.plot(tau, expectation_values[:,0], label='Sz', lw=2, alpha=.8)
plt.plot(tau, expectation_values[:,1], label='Sz', lw=2, alpha=.8)
plt.ylabel('<Mz>')
plt.xlabel('Tau [time]')
plt.xlim(tau[0], tau[-1])
plt.show()