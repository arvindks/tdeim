# Tensor Sensor Placement


This repository contains MATLAB codes for determining optimal sensor placements using tensor decompositions and accompanies the paper


"Tensor-based flow reconstruction from optimally located sensor measurements" -- M. Farazmand and A.K. Saibaba.



### Requirements
To use the package you need MATLAB 2020b or higher (the svdsketch is required for the 3D reconstructions, otherwise a lower version maybe sufficient). You will also need the [Tensor Toolbox](https://www.tensorlab.net/).

### Datasets
To run the codes, you will need to download the following datasets:
1. Sea surface temperature. [Description](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html). [data](https://psl.noaa.gov/repository/entry/show?entryid=12159560-ab82-48a1-b3e4-88ace20475cd). 
2. 3D Unsteady Research Vessel Tangaroa. [Description](https://cgl.ethz.ch/research/visualization/data.php). [data](https://cgl.ethz.ch/Downloads/Data/ScientificData/tangaroa3d_nc.zip). 

### License
To use these codes in your research, see the [license](https://github.com/arvindks/tdeim/blob/main/License.md). If you use this code in any way, please cite our paper


```
@article{farazmand2022tensor,
  title={Tensor-based flow reconstruction from optimally located sensor measurements},
  author={Farazmand, Mohammad and Saibaba, Arvind K},
  journal={arXiv preprint arXiv:2208.09875},
  year={2022}
}
```

### Funding
This work was supported, in part, by the National Science Foundation through the awards DMS-1821149 and DMS-1745654.
