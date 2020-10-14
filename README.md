

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4045184.svg)](https://doi.org/10.5281/zenodo.4045184)


# Abstract
Surface gravity waves play a major role in the exchange of momentum, heat, energy, and gases between the ocean and the atmosphere. The interaction between currents and waves can lead to variations in the wave direction, frequency, and amplitude. In the present work, we use an ensemble of synthetic currents to force the wave model WAVEWATCH III and assess the relative impact of current divergence and vorticity in modifying several properties of the waves, including direction, period, directional spreading, and significant wave height (Hs). We find that the spatial variability of Hs is highly sensitive to the nature of the underlying current and that refraction-caused vorticity in the rotational component of the flow is the main mechanism leading to gradients of Hs. The results obtained using synthetic currents were used to interpret the response of surface waves to realistic currents by running an additional set of simulations using the llc4320 MITgcm output in the California Current region. Our findings suggest that wave parameters could be used to detect and characterize strong gradients in the velocity field, which is particularly relevant for the Surface Water and Ocean Topography (SWOT) satellite as well as several proposed satellite missions.
# Authors
* [Bia Villas Boas](https://scripps.ucsd.edu/profiles/avillasboas) <<avillasboas@ucsd.edu>>
* [Sarah T. Gille](http://www-pord.ucsd.edu/~sgille/)
* [Matthew R. Mazloff](http://scrippsscholars.ucsd.edu/mmazloff)
* [Bruce D. Cornuelle](http://scrippsscholars.ucsd.edu/bcornuelle)
* [Fabrice Ardhuin](https://annuaire.ifremer.fr/cv/16811/en/)
# Data
This project uses synthetic currents created using the code. 

# WAVEWATCH III configuration

The model used in this paper was compiled with the following switches:

F90 NOGRB NOPA LRB4 SCRIP SCRIPNC NC4 TRKNC DIST MPI PR3 UQ FLX0 LN1 ST4 STAB0 NL1 BT0 DB0 TR0 BS0 IC0 IS0 REF0 IG0 XX0 WNT2 WNX1 RWND CRT1 CRX1 O0 O1 O2 O2a O2b O2c O3 O4 O5 O6 O7

# Funding
This project was partlially funded by the SWOT program with NASA grant NNX16AH67G.
Bia Villas Boas was partially funded by NASA grant 80NSSC17K0326.
