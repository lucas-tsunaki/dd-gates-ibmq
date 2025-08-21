# Multi-Qubit Gates by Dynamical Decoupling of Central Qubit Performed with IBMQ

This repository contains the python script for the graphical user interface (GUI) application of the dataset from the paper *Multi-Qubit Gates by Dynamical Decoupling of Central Qubit Performed with IBMQ and 15NV Center in Diamond* by *L. Tsunaki et al*, available at: <>

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