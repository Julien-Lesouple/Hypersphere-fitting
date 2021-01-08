# Matlab Routine
This folder containes matlab script to estimate hypersphere using EM and hel repoducing the results of the associated paper.
## main_2D
Use to estimate circle and compare with Kasa, E-Landau and ILS. Setting kappa = 0 leads to a uniform prior, and kappa>0 to a von Mises-Fisher prior
## main_3D
Use to estimate sphere and compare with FGFA and ILS. Setting kappa = 0 leads to a uniform prior, and kappa>0 to a von Mises-Fisher prior
## main_comp_time_2d
Estimates the mean computation times and associated standard deviations of each 2D algorithm
## main_comp_time_3d
Estimates the mean computation times and associated standard deviations of each 3D algorithm
## Higher dimensions
To estimate higher dimension hyperspheres please use the same main but remove variables associated with state-of-the-art algorithm and set d to whatever you want
## Remarks
vmrand was downloaded on Mathworks File exchange (https://fr.mathworks.com/matlabcentral/fileexchange/37241-vmrand-fmu-fkappa-varargin)
SphericalDistributionsRand was downloaded on Mathworks File exchange (https://fr.mathworks.com/matlabcentral/fileexchange/52398-sphericaldistributionsrand)
