# pyhesdm
*Local Universe Dispersion Measure Model Computed from HESTIA Simulation*    

### Download   
```
git clone git@github.com:yuxinhuang1229/pyhesdm.git
```      
### Test installation   
```setup.py``` automatically download NEDLVS catalog from its website. This catalog is used for calculating DM_halos  
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

# calculate DM_IGM
virgo = SkyCoord("12h27m 12d43m", frame='icrs')
hesdm.get_dmigm(virgo.galactic.l,virgo.galactic.b)
>> 29.371574 pc/cm^3

# calculate DM_halos
halos.get_dmhalos(lon=virgo.galactic.l, lat=virgo.galactic.b)
>> 618 objects are identified along the ray
>> 226 objects already have masses
>> 9  groups/clusters found within 5Mpc.
>> <DM_halos> =  71.72133087498753 pc/cm^3
>> <DM_halos>_specz =  70.33314714677843 pc/cm^3
>> <DM_halos>_photoz =  1.388183728209114 pc/cm^3
>> <DM_halos>_grp =  154.94302626129146 pc/cm^3
>> 225.27617 pc/cm^3
```

