# Non-Hermitian-Localisation

Authors: E. O. HILTUNEN and Y. DE BRUIJN

Date: 14.04.2025

------------------------------------------------------------------------------------------------------------------

In this computational notebook, we provide the MATLAB code for the computations in [1].

## I.1 Bandfunctions:

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


<p align="center"> <img src="Figures/Defected_resonator_chain.png" alt="BandMonomer" width="1000"/> </p>

<p align="center"> <img src="Figures/Schematic_Spectrum.png" alt="BandMonomer" width="600"/> </p>


<p align="center"> <img src="Figures/Defect_resonance_monomer.png" alt="BandMonomer" width="300"/> </p>

<p align="center"> <img src="Figures/Decay_lengths_monomer.png" alt="BandMonomer" width="300"/> </p>

<p align="center"> <img src="Figures/Defect_resonance_dimer.png" alt="BandMonomer" width="300"/> </p>

<p align="center"> <img src="Figures/Decay_lengths_dimer.png" alt="BandMonomer" width="300"/> </p>
