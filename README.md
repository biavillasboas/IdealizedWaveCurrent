

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4045184.svg)](https://doi.org/10.5281/zenodo.4045184)

# Source code for 
Villas Bôas, A. B., B. D. Cornuelle, M. R. Mazloff, S. T. Gille, and F. Ardhuin, Wave-Current Interactions at Meso and Submesoscales: Insights from Idealized Numerical Simulations. J. Phys. Oceanogr., doi: https://doi.org/10.1175/JPO-D-20-0151.1. 

# Abstract
Surface gravity waves play a major role in the exchange of momentum, heat, energy, and gases between the ocean and the atmosphere. The interaction between currents and waves can lead to variations in the wave direction, frequency, and amplitude. In the present work, we use an ensemble of synthetic currents to force the wave model WAVEWATCH III and assess the relative impact of current divergence and vorticity in modifying several properties of the waves, including direction, period, directional spreading, and significant wave height (Hs). We find that the spatial variability of Hs is highly sensitive to the nature of the underlying current and that refraction-caused vorticity in the rotational component of the flow is the main mechanism leading to gradients of Hs. The results obtained using synthetic currents were used to interpret the response of surface waves to realistic currents by running an additional set of simulations using the llc4320 MITgcm output in the California Current region. Our findings suggest that wave parameters could be used to detect and characterize strong gradients in the velocity field, which is particularly relevant for the Surface Water and Ocean Topography (SWOT) satellite as well as several proposed satellite missions.
# Authors
* [Bia Villas Boas](https://scripps.ucsd.edu/profiles/avillasboas) <<avillasboas@ucsd.edu>>
* [Sarah T. Gille](http://www-pord.ucsd.edu/~sgille/)
* [Matthew R. Mazloff](http://scrippsscholars.ucsd.edu/mmazloff)
* [Bruce D. Cornuelle](http://scrippsscholars.ucsd.edu/bcornuelle)
* [Fabrice Ardhuin](https://annuaire.ifremer.fr/cv/16811/en/)

# Data
All model output analyzed in this paper is availabe for download here https://doi.org/10.6075/J0X928V6

# Funding
This project was partlially funded by the [SWOT](https://swot.jpl.nasa.gov/) program with NASA grant NNX16AH67G and 80NSSC20K1136.
Bia Villas Bôas was also partially funded by NASA grant 80NSSC17K0326.

# How to use this repository

All figures in Villas Bôas et al. (2020) can be reproduced using the Python scripts from this repository and the [model output](https://library.ucsd.edu/dc/collection/bb6217292w). To do so, follow these steps

1. Make a local copy of this repository by either cloning or downloading it.

2. Download the [model output](https://library.ucsd.edu/dc/collection/bb6217292w), untar the files, and move all three directories to `data` in the project root. After doing so, your directory three should look like this:

```
IdealizedWaveCurrent/
├── data
│   ├── llc4320
│   ├── model_stats
│   └── synthetic
├── figs
├── src
└── tools
```
3. Make sure that you create an environment with the package versions specified in `environment.yml`. If you are using [Conda](https://docs.conda.io/en/latest/) you can run 

`conda env create -f environment.yml`

from the project root.

4. If you follow the steps above you should be able to reproduce all figures, by running `python figXX.py` from the `src` directory without having to adjust any paths.

* **Note on rendering matplotlib text with LaTeX:** To ensure that the math fonts in the figures matched the fonts in the paper, the code in this repository requires a [working LaTeX installation](https://matplotlib.org/3.1.0/tutorials/text/usetex.html). If you encounter any problems with this, you may change the line matplotlib.rcParams['text.usetex']=True  to False or just comment it out. 

# How to cite this code

If you wish to use the code from this repository, you may cite it as 

Villas Bôas, Ana B. (2020, September 23). Source code for: 'Wave-Current Interactions at Meso and Submesoscales: Insights from Idealized Numerical Simulations'. Zenodo. https://doi.org/10.5281/zenodo.4045183
