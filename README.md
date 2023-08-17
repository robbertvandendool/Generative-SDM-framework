# Generative-SDM-framework
This repository contains the data and code necessary to reproduce the paper van den Dool et al. (2023) 'Species distribution modelling using a generative framework' by Robbert T. van den Dool, Alejandro Morales, Wopke van der Werf & Jacob C. Douma.

# Scripts
The scripts folder contains the scripts to reproduce the results described in the: 1) performance benchmark, 2) case study 1 on Bradypus variegatus and 3) case study 2 showing parametric methods on two virtual species. No data needed to be included, since all data is freely available through R packages 'disdat' and 'dismo' or in case study 2 the data is generated in a separate script. 

# Runtime
Note that some code requires considerable time (~5-6 hours) to run. PCs with more computing cores may perform the modelling task faster through parallel computing. The number of CPU threads can be adjusted in the script(s). 