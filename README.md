# Multi-Qubit Gates by Dynamical Decoupling Implemented with IBMQ

This repository contains the basic functions for executing the DD-gate with IBMQ and simulating them with Qiskit, as presented in *Multi-Qubit Gates by Dynamical Decoupling Implemented with IBMQ and 15NV Center in Diamond by *L. Tsunaki et al*, available at: [arXiv:2508.22107](https://arxiv.org/abs/2509.22107).

The functions for generating the CPMG and XYN sequences within Qiskit SDK framework are provided at the [dd_gates.core](https://github.com/lucas-tsunaki/dd-gates-ibmq/blob/main/dd_gates/core.py) module, which can be used to obtain all the results presented in the paper.
The [examples](https://github.com/lucas-tsunaki/dd-gates-ibmq/tree/main/examples) folder provides an example on how to simulate an CPMG-10 sequence with Qiskit Aer and another on how to experimentally execute it on IBMQ hardware, as in Figure 1 (c) of the paper.

If you are looking for the simulations of the DD-gate with NV centers, please refer to the [QuaCCAToo](https://qiss-hzb.github.io/QuaCCAToo/notebooks.html) tutorials section.

## Setup

To execute the scripts in the examples section, we recommend creating a `conda` virtual environment with:
```sh
conda create --name dd-gates-env python
conda activate dd-gates-env 
```
After the virtual environment has been setup, clone the repository and run from <u>inside the local cloned folder </u>
``` sh
pip install .
```
To check the installation, you can run
``` sh
pip show dd_gates
```
This application has dependencies on NumPy, Matplotlib, Qiskit version 2.1.0, Qiskit Aer Simulator version 0.17.1 and Qiskit IBM Runtime version 0.40.1.
Due to Qiskit's fast development, the code is not guaranteed to work with newer versions as the ones specified above. 
