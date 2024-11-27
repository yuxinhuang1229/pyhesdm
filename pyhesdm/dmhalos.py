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

class NEDLVS_Tully_Halos:
    
    def __init__(self, ned_file = f'{os.path.dirname(__file__)}/NEDLVS_20210922_v2.fits',
                tully_file = f'{os.path.dirname(__file__)}/tully_groups.csv'):
        self.cosmo = FlatLambdaCDM(H0=67.77, Om0=0.270, Tcmb0=2.725, Ob0=0.048)
        self.nedlvs_tab = self.get_nedlvs_tab(ned_file)
        self.tully_tab = self.get_tully_tab(tully_file)
        self.mpch2z = self.dist2z_interp(cosmology=self.cosmo, num_points=100000)
        
    def dist2z_interp(self, cosmology, num_points): # Mpc/h to z
        """
        Returns:
        Interpolation function: comoving distance (Mpc/h) -> redshift
        """
        # Create a grid of redshift values
        z_values = np.linspace(0, 0.1, num_points)  # Adjust the upper limit based on your needs

        # Calculate the corresponding comoving distances for the given redshift grid
        comoving_distances = cosmology.comoving_distance(z_values)

        # Interpolation function: comoving distance (Mpc/h) -> redshift
        interp_function = interp1d(comoving_distances.to(u.Mpc/cosmology.h).value, z_values, kind='linear', fill_value="extrapolate")

        return interp_function
         
    def get_nedlvs_tab(self, ned_file):
        """
        Returns:
        nedlvs_tab (Table): NED-LVS catalog with redshift distance re-calculated using pyhesdm cosmology.
        """
        nedlvs_tab = Table.read(ned_file)
        nedlvs_tab['coord'] = SkyCoord(nedlvs_tab['ra'], nedlvs_tab['dec'], unit="deg")
        # Set redshift distances using my cosmology
        redshift_dist_sources = (nedlvs_tab['DistMpc_method']=='Redshift')
        nedlvs_tab['DistMpc'][redshift_dist_sources] = self.cosmo.luminosity_distance(nedlvs_tab['z'][redshift_dist_sources]).to('Mpc').value
        return nedlvs_tab
        
    def get_tully_tab(self, tully_file):
        """
        Returns:
        tully_tab (Table): Tully Cluster catalog.
        """
        tully_clusters = Table.from_pandas(pd.read_csv(tully_file))
        tully_clusters_coord = SkyCoord(tully_clusters['GLong']*u.deg, tully_clusters['GLat']*u.deg, frame="galactic")
        tully_clusters['ra'] = tully_clusters_coord.icrs.ra.deg
        tully_clusters['dec'] = tully_clusters_coord.icrs.dec.deg
        tully_clusters['coord'] = SkyCoord(tully_clusters['ra'], tully_clusters['dec'],unit="deg")
        tully_clusters = tully_clusters[(tully_clusters['Ng']>4)&(tully_clusters['L_Mass12']>0)]
        return tully_clusters
    
    def halo_dm(self, z:float, offset:u.Quantity, log_mhalo:float,
                rmax:float = 1.0, fhot:Callable = lambda log_mhalo: 0.25 if log_mhalo < 13.5 else 0.8, 
                log_mass_cut=13.5, halo_model:ModifiedNFW=ModifiedNFW)->u.Quantity:
        """
        For a given halo compute the DM contribution.
        Args:
            z (float): redshift
            offset (Quantity): Transverse physical offset from FRB.
            log_mhalo (float): log halo mass
            rmax (float, optional): Maximum halo gas radius in units of virial radii.
            fhot (function, optional): Ionized fraction of halo gas as a function of log halo mass.
            log_mass_cut (float, optional): Mass threshold of galaxy and cluster Mhalo. 
            halo_model (ModifiedNFW, optional): Halo gas occupation model. ModifiedNFW or any of
                its child classes defined in frb.halos.models.
        Returns:
            DM (float): Halo DM contribution in pc/cm^3. Includes the 1+z factor.
        """
        # Is it a galaxy/group halo?
        if halo_model==ModifiedNFW:
            # Just to speed things up a bit. Don't expect
            # Galaxy groups to be larger than 1 Mpc 
            if offset >= 1*u.Mpc:
                return 0.
            else:
                kwargs = {'y0':2., 'alpha':2., 'f_hot':fhot(log_mhalo)}
        elif halo_model==ICM:
            assert log_mhalo > log_mass_cut, f"You've requested the ICM model for a halo that's less massive than 10^{log_mass_cut} Msun."
            kwargs = {'f_hot':fhot(log_mhalo)}

        # Instantiate halo
        mnfw = halo_model(log_Mhalo = log_mhalo,  z = z, **kwargs)

        # Compute halo DM
        DM = mnfw.Ne_Rperp(offset, rmax=rmax).to('pc/cm**3').value/(1+z)

        # Return
        return DM
    
    def get_dmhalos(self, lon:u.Quantity, lat:u.Quantity, source_dist:u.Quantity=100*u.Mpc,
                    galaxy_model:ModifiedNFW = ModifiedNFW, cluster_model:ModifiedNFW = ICM,
                    log_mass_cut:float = 13.5, rmax:float = 1., 
                    fhot:Callable = lambda log_mhalo: 0.25 if log_mhalo < 13.5 else 0.8, 
                    if_print:bool=True, full_output:bool=False)->tuple:
        """
        Estimate the contribution of the intersecting halos listed
        within the NEDLVS catalog.
        Warning: Ignores halos without mass estimates in the table.

        Args:
            lon (u.Quantity): Longitude of the source
            lat (u.Quantity): Latitude of the source
            source_dist: Distance from the source to the observer, 3.4Mpc < source_dist < 120Mpc
            galaxy_model (ModifiedNFW, optional): Halo gas occupation model of galaxies. ModifiedNFW or any of
                its child classes defined in frb.halos.models.
            cluster_model (ModifiedNFW, optional): Halo gas occupation model of clusters. ICM or any of
                its child classes defined in frb.halos.models.
            log_mass_cut (float, optional): Mass threshold of galaxy and cluster Mhalo.
            rmax (float, optional): Maximum halo gas radius in units of virial radii.
            fhot (function, optional): Ionized fraction of halo gas as a function of log halo mass.
            if_print (bool, optional): If True, print halo information along the sightline.
            full_output (bool, optional): If True, output all components of DM_halos and list galaxies and groups; 
                else only print the final DM_halos.
        Returns:
            if full_output is True:
                mean_dm_halos_lvs (float): Mean value of DM_halos from galaxies with all the methods.
                mean_dm_halos_lvs_spec (float): Mean value of DM_halos from galaxies with spec-zs 
                                                and redshift-independent distances.
                mean_dm_halos_lvs_phot (float): Mean value of DM_halos from galaxies with photo-z's and other methods.
                mean_grp_dm (float): Mean value of DM_halos from the Tully cluster catalog.
                close_by_withmass (Table): Table of galaxy halos within 1 Mpc having masses from the NED LVS.
                match_grps (Table): Table of intersecting Tully Cluster catalog members. 
            if full_output is False:
                mean_dm_tot: Total DM halos (mean_dm_halos_lvs_spec+mean_grp_dm)
        """
        
        if (source_dist <= 3.4*u.Mpc) or (source_dist > 120*u.Mpc):
            raise ValueError("Distance must be between 3.4Mpc and 120Mpc!")
        
        nedlvs_tab = self.nedlvs_tab.copy()
        tully_tab = self.tully_tab.copy()    
        
        ray_coord = SkyCoord(l=lon, b=lat, frame='galactic')
        #import pdb; pdb.set_trace()
        # Estimate separations of objects within the table
        nedlvs_tab['ang_sep'] = ray_coord.separation(nedlvs_tab['coord']).to('arcmin')
        nedlvs_tab['phys_sep'] = nedlvs_tab['DistMpc']*u.Mpc*np.sin(nedlvs_tab['ang_sep'].to('rad').value)

        # Conditions for closeness
        #import pdb; pdb.set_trace() 
        valid_distances = nedlvs_tab['DistMpc']>0 # Weed out weird ones
        distance_cut = nedlvs_tab['DistMpc']<=source_dist.to('Mpc').value #Only need foreground objects
        local_grp_cut = nedlvs_tab['DistMpc']>3.4 # Exclude local group
        phys_sep_cut = nedlvs_tab['phys_sep']<1*u.Mpc # Impact param within 1 Mpc
        ang_sep_cut = nedlvs_tab['ang_sep']<90*u.deg # Make sure the earth is not between the FRB and the galaxy
        z_secure = ~nedlvs_tab['z_qual'] # Exclude redshift with bad quality
        is_nearby_fg = valid_distances & distance_cut & local_grp_cut & phys_sep_cut & ang_sep_cut & z_secure

        # Get culled table
        close_by = nedlvs_tab[is_nearby_fg]
        

        # Have mass estimates?
        mass_cut = (close_by['Mstar']>0) & (close_by['Mstar']<1e13) # Weird to have galaxies above this value
        close_by_withmass = close_by[mass_cut]
        

        if len(close_by_withmass)==0:
            if if_print:
                print(f"No local volume halos with mass found along this sightline (l={lon}, b={lat}, distance={source_dist}).")
            if full_output:
                return 0.*u.pc/(u.cm**3), 0.*u.pc/(u.cm**3), 0.*u.pc/(u.cm**3), 0.*u.pc/(u.cm**3), Table(), Table()
            else:
                return 0.*u.pc/(u.cm**3)
            
        # Get DM estimates
        # Start with Mstar. Reduce by 0.3 dex for Salpeter to Chabrier IMF
        close_by_withmass['log_mstar'] = np.log10(close_by_withmass["Mstar"]/1.7)
        logmhalo_column = MaskedColumn(name='log_mhalo', dtype=float, length=len(close_by_withmass))
        close_by_withmass.add_column(logmhalo_column)
        close_by_withmass = close_by_withmass.filled(-99.)
        
        halo_dms = []
        for num, (offset, distance) in enumerate(zip(close_by_withmass['phys_sep'], close_by_withmass['DistMpc'])):
            try:
                if (close_by_withmass[num]['DistMpc_method'] == 'zIndependent') & (close_by_withmass[num]['z_tech']!='SPEC'):
                    close_by_withmass[num]['z'] = self.mpch2z(distance*self.cosmo.h)
                close_by_withmass[num]['log_mhalo'] = halomass_from_stellarmass(close_by_withmass[num]['log_mstar'], 
                                                                                z=close_by_withmass[num]['z'])
                hdm = self.halo_dm(z=close_by_withmass[num]['z'], 
                                   offset=offset*u.Mpc, 
                                   log_mhalo=close_by_withmass[num]['log_mhalo'], 
                                   rmax=rmax, fhot=fhot, 
                                   halo_model=galaxy_model)
                halo_dms.append(hdm)
            except:
                print(f"Something went wrong for halo # {num}")
                import pdb; pdb.set_trace()
        close_by_withmass['DM_halo'] = halo_dms
        close_by_withmass.sort("DM_halo",reverse=True)

        # Are there any group/cluster members in the field?
        # Compute transverse distances and initialize columns
        tully_cluster_centers = tully_tab[np.isin(tully_tab['pgc'], np.unique(tully_tab['1PGC']))]
        match_grps = tully_cluster_centers.copy()
        match_grps['ang_sep'] = ray_coord.separation(match_grps['coord']).to('arcmin')
        match_grps['phys_sep'] = match_grps['D_v'].value/self.cosmo.h*u.Mpc*np.sin(match_grps['ang_sep'].to('rad').value)
        grps_infield = (match_grps['phys_sep']<5*u.Mpc)&(match_grps['ang_sep']<90*u.deg)&(match_grps['D_v'].value/self.cosmo.h<=source_dist.to('Mpc').value)
        match_grps = match_grps[grps_infield]

        # Cross match with the tully cluster catalog
        if len(match_grps)>0:
            use_grps = True
            match_grps['DM_halo'] = 0.0
            match_lvs, _ = xmatch_catalogs(close_by_withmass, tully_tab, skydist=1*u.arcsec)

            if len(match_lvs)>0:
                # remove matched objects from the table
                close_by_withmass = setdiff(close_by_withmass, match_lvs, keys="objname")
        else:
            use_grps = False

        
        # Loop through groups
        if use_grps:
            for central_entry in match_grps:
                # Uses luminosity-based mass but change this to virial if needed.
                log_grp_mass = np.log10(central_entry['L_Mass12'])+12.

                # Use the ICM model if log halo mass is above log_mass_cut
                if log_grp_mass <= log_mass_cut:
                    grp_model = galaxy_model
                else:
                    grp_model = cluster_model
                central_entry['DM_halo'] = self.halo_dm(z=self.mpch2z(central_entry['D_v']),
                                      offset=central_entry['phys_sep']*u.Mpc,
                                      log_mhalo=log_grp_mass, rmax=rmax, fhot=fhot,
                                      halo_model=grp_model)
        else:
            match_grps = Table()


        # return average dm_halos
        mean_dm_halos_lvs = np.sum(close_by_withmass["DM_halo"])
        if use_grps:
            mean_grp_dm = np.sum(match_grps['DM_halo'])
        else:
            mean_grp_dm = 0.0
        spec_z_halos = (close_by_withmass['z_tech']=='SPEC')|(close_by_withmass['DistMpc_method']=='zIndependent')
        mean_dm_halos_lvs_spec = np.sum(close_by_withmass['DM_halo'][spec_z_halos])
        mean_dm_halos_lvs_phot = np.sum(close_by_withmass['DM_halo'][~spec_z_halos])
        
        if if_print:
            print(f"{len(close_by)} objects are identified along the ray")
            print(f"{len(close_by_withmass)} objects already have masses")
            print(len(match_grps), " groups/clusters found within 5Mpc.")
            print("<DM_halos> = ", mean_dm_halos_lvs, "pc/cm^3")
            print("<DM_halos>_specz = ", mean_dm_halos_lvs_spec, "pc/cm^3")
            print("<DM_halos>_photoz = ", mean_dm_halos_lvs_phot, "pc/cm^3")
            print("<DM_halos>_grp = ", mean_grp_dm, "pc/cm^3")
        
        if full_output:
            return mean_dm_halos_lvs*u.pc/(u.cm**3), mean_dm_halos_lvs_spec*u.pc/(u.cm**3), mean_dm_halos_lvs_phot*u.pc/(u.cm**3), mean_grp_dm*u.pc/(u.cm**3), close_by_withmass, match_grps
        else:
            mean_dm_tot = (mean_dm_halos_lvs_spec+mean_grp_dm)*u.pc/(u.cm**3)
            return mean_dm_tot    
        
        
        