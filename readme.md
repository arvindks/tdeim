# Tensor Sensor Placement


This repository contains MATLAB codes for determining optimal sensor placements using tensor decompositions and accompanies the paper

Farazmand, Mohammad, and Arvind K. Saibaba. "Tensor-based flow reconstruction from optimally located sensor measurements." Journal of Fluid Mechanics 962 (2023): A27. [Website](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/tensorbased-flow-reconstruction-from-optimally-located-sensor-measurements/E6F6DD6727238BC50C42919936E4D841). [arxiv](https://arxiv.org/abs/2208.09875).



### Requirements
To use the package you need MATLAB 2020b or higher (the svdsketch is required for the 3D reconstructions, otherwise a lower version maybe sufficient). You will also need the [Tensor Toolbox](https://www.tensorlab.net/).

### Datasets
To run the codes, you will need to download the following datasets:
1. Kolmogorov flow. [data](https://doi.org/10.5281/zenodo.7464956).
2. Sea surface temperature. [Description](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html). [data](https://psl.noaa.gov/repository/entry/show?entryid=12159560-ab82-48a1-b3e4-88ace20475cd). 
3. 3D Unsteady Research Vessel Tangaroa. [Description](https://cgl.ethz.ch/research/visualization/data.php). [data](https://cgl.ethz.ch/Downloads/Data/ScientificData/tangaroa3d_nc.zip). 

### License
To use these codes in your research, see the [license](https://github.com/arvindks/tdeim/blob/main/License.md). If you use this code in any way, please cite our paper


```
@article{farazmand_saibaba_2023, title={Tensor-based flow reconstruction from optimally located sensor measurements}, volume={962}, DOI={10.1017/jfm.2023.269}, journal={Journal of Fluid Mechanics}, publisher={Cambridge University Press}, author={Farazmand, Mohammad and Saibaba, Arvind K.}, year={2023}, pages={A27}}
```

### Funding
This work was supported, in part, by the National Science Foundation through the awards DMS-1821149 and DMS-1745654.
