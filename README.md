# Spin Half Quantum Circuit Program

This program is designed to simulate the dynamics (including entanglement entropy and mutual information) of a quantum circuit of spin-1/2 particles
evolving subject to local disorder and continuous measurements. 

Dr Marcin Szyniszewski and I developed this program as part of our publication,
Diagnostics of entanglement dynamics in noisy and disordered spin chains via the measurement-induced steady-state. 

Phys. Rev. B 105, 144202.

## Usage

The program takes shell inputs (described below) and returns up to three outputs,
- Entanglement entropy
- Bi-partite mutual information (if system length is divisible by 4)
- Tri-partite mutual information (if system length is divisible by 4)

### Inputs for shell script:

- circuitLength   (Number of qubits)
- piece           (See: How mutual information (MI) works)
- measProbability (Probability for measurement, program scales this with dt)
- measStrength    (Measurement Strength)
- xiStatic        (Strength of static disorder)
- xiTemporal      (Strength of uniform temporal disorder, program scales this with dt)
- xiNonstatic     (Strength of non-static random disorder, program scales this with dt)
- dt              (Time step increment)
- startReadings   (When to start entropy readings - real time)
- entropyTiming   (Timing between entropy readings in units of dt) older versions used units of real time but this led to problematic floating point errors
- simTime         (Time over which sim is performed - real time)

Description of mutual information calculation can be found in the preamble of the program.  

## Requirements
- Eigen package
