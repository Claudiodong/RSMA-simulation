# RSMA Simulation and Wireless Communication Experiments

This repository contains MATLAB simulations for wireless communication and optimisation experiments. 

## Project Structure

| Folder | Contents |
| --- | --- |
| `practice/` | Water-filling practice implementation used to explore power allocation across parallel channels. |
| `cw1-mimo-capacity-and-ber/` | MIMO capacity, water-filling, QPSK link simulations, and BER comparisons across SISO, SIMO, MISO, and MIMO configurations. |
| `cw2-massive-mimo-simulation/` | Massive MIMO simulations covering random user deployment, long-term SINR, channel time correlation, scheduling, zero-forcing beamforming, and block diagonalisation experiments. |
| `cw3-beamforming-and-ris-optimisation/` | Beamforming and RIS-aided optimisation experiments using KKT, SDP, SOCP, WMMSE, and zero-forcing/water-filling methods. |

## Result Gallery

### MIMO Capacity and BER

<p align="center">
  <img src="cw1-mimo-capacity-and-ber/results/Capacity.jpg" width="360" alt="CW1 capacity result">
</p>

<p align="center">
  <img src="cw1-mimo-capacity-and-ber/results/ALL%20results.jpg" width="360" alt="CW1 BER comparison">
</p>

<p align="center">
  <img src="cw1-mimo-capacity-and-ber/results/MIMO%20SM%20reception.jpg" width="360" alt="CW1 MIMO spatial multiplexing">
</p>

### Massive MIMO Simulation

<p align="center">
  <img src="cw2-massive-mimo-simulation/results/Random_deployment.png" width="360" alt="CW2 random deployment">
</p>

<p align="center">
  <img src="cw2-massive-mimo-simulation/results/Long_term%20SINR.png" width="360" alt="CW2 long-term SINR">
</p>

<p align="center">
  <img src="cw2-massive-mimo-simulation/results/Task6_result.png" width="360" alt="CW2 block diagonalisation result">
</p>

### Beamforming and RIS Optimisation

<p align="center">
  <img src="cw3-beamforming-and-ris-optimisation/results/dB_optimisation_plot_task1.png" width="360" alt="CW3 beamforming optimisation">
</p>

<p align="center">
  <img src="cw3-beamforming-and-ris-optimisation/results/Task2.png" width="360" alt="CW3 WMMSE comparison">
</p>

<p align="center">
  <img src="cw3-beamforming-and-ris-optimisation/results/Task4_increase_M.png" width="360" alt="CW3 RIS antenna sweep">
</p>

## Running the Simulations

Open MATLAB from the repository root and add the relevant source folders before running an experiment script. For example:

```matlab
addpath(genpath("cw1-mimo-capacity-and-ber"));
run("cw1-mimo-capacity-and-ber/main.m");
```

For the massive MIMO and optimisation experiments, use the corresponding top-level scripts inside `cw2-massive-mimo-simulation/legacy`, `cw2-massive-mimo-simulation/latest`, or `cw3-beamforming-and-ris-optimisation/`.

Author: Claudio Dong
