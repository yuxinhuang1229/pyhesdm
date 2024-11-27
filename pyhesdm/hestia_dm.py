import numpy as np
import pandas as pd
import healpy as hp
import os

from astropy import units as u
from astropy.io import fits
from astropy.table import Table, vstack, join, setdiff, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.coordinates import Longitude
from astropy.cosmology import FlatLambdaCDM,Cosmology,z_at_value
from astropy.stats import sigma_clipped_stats

from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.optimize import fsolve

from typing import Union, Callable

from frb.halos.models import ModifiedNFW, halomass_from_stellarmass, ICM, YF17
from frb import mw
from frb.dm.igm import average_DMIGM
from frb.frb import FRB
from frb.galaxies.frbgalaxy import FRBHost
from frb.surveys.catalog_utils import xmatch_catalogs

class Hestia_DM:
    
    "DM model of Local Universe Calculated from HESTIA Simulation"
    
    def __init__(self, model_csv=f'{os.path.dirname(__file__)}/model_lg_nside64.csv', nside=64):
        self.nside = nside
        self.model = pd.read_csv(model_csv)
        self.dwarfs = pd.read_csv(f'{os.path.dirname(__file__)}/dwarfs.csv')
        self.igm_model = pd.read_csv(f'{os.path.dirname(__file__)}/dmigm_layers_8Mpc.csv')
        self.igm_std = pd.read_csv(f'{os.path.dirname(__file__)}/dmigm_std_layers_8Mpc.csv')
    
    def galactic_to_healpix(self, lon:u.Quantity, lat:u.Quantity):
        """
        Convert galactic latitude and longitude to HEALPix pixel number.
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
        Returns:
            pixel (int or array of int): HEALPix pixel(s) corresponding to the given Galactic Longitude and Latitude.
        """
        theta = np.pi/2 - lat.to(u.radian).value  # convert latitude to colatitude
        phi = lon.to(u.radian).value         # convert longitude to radians
        pixel = hp.ang2pix(self.nside, theta, phi, nest=False)
        return pixel
    
    def get_disk(self, lon:u.Quantity, lat:u.Quantity):
        """
        Compute the disk DM from NE2001 model.
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
        Returns:
            dm_ne2001 (astropy.units.Quantity): DM of the NE2001 component.
        """
        pixel = self.galactic_to_healpix(lon, lat)
        dm_ne2001 = np.array(self.model.loc[pixel,'ne2001']) * u.pc/(u.cm**3)
        return dm_ne2001
    
    def get_mwhalo(self, lon:u.Quantity, lat:u.Quantity):
        """
        Compute the Milky Way halo DM from HESTIA simulation.
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
        Returns:
            dm_mwhalo (astropy.units.Quantity): DM of the Milky Way halo component.
        """
        pixel = self.galactic_to_healpix(lon, lat)
        self.in_dwarf_halo(lon, lat, dwarfs=[0,1])
        dm_mwhalo = np.array(self.model.loc[pixel,'mw_halo']) * u.pc/(u.cm**3)
        return dm_mwhalo
    
    def get_lgigrm(self, lon:u.Quantity, lat:u.Quantity):
        """
        Compute the Local Group intra-group medium DM from HESTIA simulation.
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
        Returns:
            dm_lgigrm (astropy.units.Quantity): DM of the Local Group Intra-group Medium component.
        """
        pixel = self.galactic_to_healpix(lon, lat)
        self.in_dwarf_halo(lon, lat, dwarfs=[2,3])
        dm_lgigrm = np.array(self.model.loc[pixel,'lg_igrm']) * u.pc/(u.cm**3)
        return dm_lgigrm
    
    def get_mw(self, lon:u.Quantity, lat:u.Quantity):
        """
        Compute the total Milky Way DM (DM_halo+DM_disk).
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
        Returns:
            dm_mw (astropy.units.Quantity): DM of the Milky Way (including NE2001 and halo) component.
        """
        pixel = self.galactic_to_healpix(lon, lat)
        self.in_dwarf_halo(lon, lat, dwarfs=[0,1])
        dm_mw = (np.array(self.model.loc[pixel,'mw_halo']) 
                + np.array(self.model.loc[pixel,'ne2001'])) * u.pc/(u.cm**3)
        return dm_mw
    
    def get_lg(self, lon:u.Quantity, lat:u.Quantity):
        """
        Compute the total Local Group DM (DM_halo+DM_disk+DM_igrm).
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
        Returns:
            dm_lg (astropy.units.Quantity): DM of the Local Group (including NE2001, halo and Intra-group Medium) component.
        """
        pixel = self.galactic_to_healpix(lon, lat)
        self.in_dwarf_halo(lon, lat)
        dm_lg = (np.array(self.model.loc[pixel,'lg_igrm']) + np.array(self.model.loc[pixel,'mw_halo']) 
                + np.array(self.model.loc[pixel,'ne2001'])) * u.pc/(u.cm**3)
        return dm_lg
    
    def in_dwarf_halo(self, lon:u.Quantity, lat:u.Quantity, dwarfs:list=[0,1,2,3]):
        """
        Whether the sightline is going through any dwarf galaxy halo (LMC, SMC, M31, M33)
        """
        targ_num = 1
        if lon.isscalar:
            lons = [lon]
            lats = [lat]
        else:
            lons = lon
            lats = lat
        for l,b in zip(lons,lats):
            skycoord = SkyCoord(l,b,frame='galactic')
            dwarfs_coord = SkyCoord(self.dwarfs.loc[dwarfs,'Lon'],self.dwarfs.loc[dwarfs,'Lat'], unit='deg', frame='galactic')
            self.dwarfs.loc[dwarfs,'separation'] = skycoord.separation(dwarfs_coord).to("deg").value
            self.dwarfs['in_dwarf_halo']=False
            self.dwarfs.loc[dwarfs,'in_dwarf_halo'] = (self.dwarfs.loc[dwarfs,'separation']<self.dwarfs.loc[dwarfs,'Angle'])
            if any(self.dwarfs.loc[dwarfs,'in_dwarf_halo']):
                print(f'Target Number {targ_num}: Longitude={l},Latitude={b}')
                for name in self.dwarfs.loc[self.dwarfs['in_dwarf_halo'],'Name']:
                    print(f'This sightline passes the {name} halo!')
            targ_num += 1
            
    def get_dmigm_1line(self, lon:u.Quantity, lat:u.Quantity, source_dist:u.Quantity=100*u.Mpc, figm=0.8):
        """
        Compute the DM_igm for one sightline
        source_dist should be within 3.4Mpc < source_dist < 120Mpc
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can only be a scalar.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can only be a scalar.
            source_dist (astropy.units.Quantity, optional): Comoving distance of the source, provided as an Astropy Quantity. \
            Can only be a scalar.
            figm (float, optional): the fraction of cosmic baryons in IGM
        Returns:
            DM of the IGM component and its standard deviation calculated from all the realizations.
        """
        pixel = self.galactic_to_healpix(lon, lat)
        dist_Mpc = source_dist.to(u.Mpc).value
        dist_samples = np.hstack((np.arange(4,120,8),np.array([120])))
        dm_samples = self.igm_model.iloc[pixel].to_numpy() * figm/0.8
        dm_std_samples = self.igm_std.iloc[pixel].to_numpy() * figm/0.8
        return np.interp(dist_Mpc, dist_samples, dm_samples), np.interp(dist_Mpc, dist_samples, dm_std_samples)
    
    def get_dmigm(self, lon:u.Quantity, lat:u.Quantity, source_dist:u.Quantity=100*u.Mpc, figm=0.8):
        """
        Compute the DM_igm for multiple sightlines
        source_dist should be within 3.4Mpc < source_dist < 120Mpc
        Args: 
            lon (astropy.units.Quantity): Galactic Longitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            lat (astropy.units.Quantity): Galactic Latitude of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            source_dist (astropy.units.Quantity, optional): Comoving distance of the source, provided as an Astropy Quantity. \
            Can be a scalar or an array.
            figm (float, optional): the fraction of cosmic baryons in IGM
        Returns:
            DM of the IGM component and its standard deviation calculated from all the realizations.
        """
        if (source_dist <= 3.4*u.Mpc) or (source_dist > 120*u.Mpc):
            raise ValueError("Distance must be between 3.4Mpc and 120Mpc!")
        dmigm_list = []
        dmstd_list = []
        if lon.isscalar:
            dmigm, dmstd = self.get_dmigm_1line(lon,lat,source_dist,figm=figm)
            return dmigm*u.pc/(u.cm**3), dmstd*u.pc/(u.cm**3)
        else:
            lons = lon
            lats = lat
            if type(source_dist)!=list:
                dist_Mpc = source_dist.to(u.Mpc).value
                dists = [dist_Mpc]*len(lons)*u.Mpc
            else:
                dists = source_dist
            for l,b,dist in zip(lons,lats,dists):
                dmigm, dmstd = self.get_dmigm_1line(l,b,dist,figm=figm)
                dmigm_list.append(dmigm)
                dmstd_list.append(dmstd)
            return np.array(dmigm_list) * u.pc/(u.cm**3), np.array(dmstd_list) * u.pc/(u.cm**3)
            
        
        