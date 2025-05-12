# Complex Band Structure for tridiagonal k-Toeplitz operators and subwavelength localisation in defected finite non-Hermitian systems

Authors: E. O. HILTUNEN and Y. DE BRUIJN

Date: 15.05.2025

------------------------------------------------------------------------------------------------------------------

In this computational notebook, we provide the MATLAB code for the computations in [1].

## Complex Quasiperiodic Gauge Capacitance matrix:

### Subwavelength Regime
In a setting where $0 < \delta \ll 1$, we seek non-trivial eigenfrequencies for the Helmholtz scattering problem such that $\omega(\delta)\to 0$, as $\delta \to 0$.

#### Quasiperiodic Capacitance
In the subwavelength regime, the complex quasiperiodic capacitance matrix enables us to derive explicit formulas for the band and gap functions.
- `Monomer_Band_Surface.m`

<p align="center"> <img src="Figures/Band_surface.png" alt="BandMonomer" width="600"/> </p>


- `Monomer_Band.m`
<p align="center"> <img src="Figures/Band_Monomer.png" alt="BandMonomer" width="300"/> </p>

- `Dimer_Band.m`
<p align="center"> <img src="Figures/Band_Dimer.png" alt="BandMonomer" width="300"/> </p>

- `Regions_monomer_Band.m`
<p align="center"> <img src="Figures/Schematic_Spectrum.png" alt="BandMonomer" width="600"/> </p>


## II. Defected Resonator chians  
A defected monomer chain of $N$ one-dimensional subwavelength resonators, with length $\ell$ and spacing $s$. The wave speed inside the resonators is $v = 1$ whereas in the defected resonator the speed is $\tilde{v} = 1 + \eta.$
<p align="center"> <img src="Figures/Defected_resonator_chain.png" alt="BandMonomer" width="1000"/> </p>

- `Defected_Spectral_Convergence.m`
<p align="center"> <img src="Figures/Defect_resonance_monomer.png" alt="BandMonomer" width="300"/> </p>

- `Monomer_Band.m`
<p align="center"> <img src="Figures/Decay_lengths_monomer.png" alt="BandMonomer" width="300"/> </p>

- `Defected_Dimers.m`
<p align="center"> <img src="Figures/Defect_resonance_dimer.png" alt="BandMonomer" width="300"/> </p>

- `Dimer_Band.m`
<p align="center"> <img src="Figures/Decay_lengths_dimer.png" alt="BandMonomer" width="300"/> </p>

## III. Non-Hermitian tight-binding Hamiltonian 

- `Defect_modes_Non_Hermitian_Hamiltonian.m`
<p align="center"> <img src="Figures/Regions_Non_Hermitian_Hamiltonian.png" alt="BandMonomer" width="300"/> </p>


## IV. References:
When using the code in the repository, please cite the following reference:

[1] De Bruijn, Y. and Hiltunen, E.O. (2025), *Complex Band Structure for tridiagonal k-Toeplitz operators and subwavelength localisation in defected finite non-Hermitian systems*.
