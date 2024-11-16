# pyhesdm
*Local Universe Dispersion Measure (DM) Model Computed from HESTIA Simulation*   

This package provides a Python interface to a DM model for the local universe within 120 Mpc, as introduced in [Huang et al.](http://arxiv.org/abs/2410.22098)(2024).    

This package is designed to calculate the DM contribution from specific components for given Galactic coordinates. The model considers the following components: the Milky Way disk, the Milky Way halo, the intra-group medium (IGrM), the intergalactic medium (IGM), and the halos. The DM models for the Milky Way disk, Milky Way halo (D < 0.2 Mpc) and the IGrM (0.2 Mpc < D < 3.4 Mpc) are all-sky projected model, while the DM model for the IGM and halos (3.4 Mpc < D < 120 Mpc) are 3D spherical models.    

The Milky Way disk DM model is the NE2001 model. If you require the NE2001 model in 3D and with higher resolution, please use [pygedm](https://github.com/FRBs/pygedm) (doi:[10.1017/pasa.2021.33](https://ui.adsabs.harvard.edu/abs/2021PASA...38...38P/abstract)). The DM models for the Milky Way halo and the IGrM are computed from the HESTIA simulation ([Libeskind et al.](https://ui.adsabs.harvard.edu/abs/2020MNRAS.498.2968L), 2020). The DM model for the IGM is computed based on the Hamlet reconstructed density field ([Valade et al.](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5148V), 2022). The DM model for the galaxy and group halos are from the NEDLVS catalog (doi:[10.26132/NED8](https://catcopy.ipac.caltech.edu/dois/doi.php?id=10.26132/NED8)) and the 2MASS Galaxy Groups catalog ([Tully et al.](https://ui.adsabs.harvard.edu/abs/2015AJ....149..171T), 2015).   

### Requirements

Before using this package, ensure that the following packages are installed:

1. [ne2001-0.0.1](https://pypi.org/project/ne2001/)   
2. [dust_extinction-1.5](https://pypi.org/project/dust-extinction/)   
3. [linetools-0.3.2](https://pypi.org/project/linetools/)
4. [mwprop-1.0.10](https://github.com/stella-ocker/mwprop)
5. [FRB](https://github.com/FRBs/FRB)   
6. [importlib_resources](https://pypi.org/project/importlib-resources/)
7. [healpy-1.18.0](https://pypi.org/project/healpy/)

### Download   
```
git clone git@github.com:yuxinhuang1229/pyhesdm.git
```      
### Test installation   
```setup.py``` automatically download NEDLVS catalog from its website. This catalog is used for calculating DM_halos.  
```
python3 setup.py sdist bdist_wheel
pip3 install .
```   
### Usage
```python
from pyhesdm.hestia_dm import Hestia_DM
from pyhesdm.hestia_dm import NEDLVS_Tully_Halos

hesdm = Hestia_DM()
halos = NEDLVS_Tully_Halos()

# calculate MW halo DM
FRB20220319D = SkyCoord('02:08:42.7 +71:02:06.9', unit=['hourangle','deg'], frame='icrs')
hesdm.get_mwhalo(FRB20220319D.galactic.l,FRB20220319D.galactic.b)
>> 39.368945 pc/cm^3

# calculate DM_IGM and its uncertainty
virgo = SkyCoord("12h27m 12d43m", frame='icrs')
dmigm, dmstd = hesdm.get_dmigm(virgo.galactic.l,virgo.galactic.b)
print(dmigm, dmstd)
>> 24.231488 pc/cm^3  2.021273 pc/cm^3

# calculate DM_halos
halos.get_dmhalos(lon=virgo.galactic.l, lat=virgo.galactic.b)
>> 618 objects are identified along the ray
>> 226 objects already have masses
>> 4  groups/clusters found within 5Mpc.
>> <DM_halos> =  71.72133087498753 pc/cm^3
>> <DM_halos>_specz =  70.33314714677843 pc/cm^3
>> <DM_halos>_photoz =  1.388183728209114 pc/cm^3
>> <DM_halos>_grp =  55.07165041895517 pc/cm^3
>> 125.4048 pc/cm^3
```

