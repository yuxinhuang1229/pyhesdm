# pyhesdm
Local Universe Dispersion Measure Model Computed from HESTIA Simulation    
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
```
