import yt
from yt.mods import *
import numpy as np
from scipy.interpolate import interp1d
#
# 1D spherical radius
#
def _sph_radius(field,data):
    return data['r']

def _sph_volume(field,data):
    #return 4.0*np.pi*data["dr"]*(data["r"]**2)
    rout = data["r"]+0.5*data["dr"]
    rin  = data["r"]-0.5*data["dr"]
    return 4.0/3.0*np.pi*(rout**3-rin**3)

def _sph_cell_mass(field,data):
    return data["dens"]*data["sph_cell_volume"]

def _sph_radial_velocity(field,data):
    return data["velx"]
def _sph_tangential_velocity(field,data):
    return data["velx"]*0.0

#
# 2D cylindrical radius
#

def _cyl_radius(field,data):
    return np.sqrt(data['r']**2 + data['z']**2)

def _cyl_volume(field,data):
    return data["dr"]*data["dz"]*2*np.pi*data["r"]

def _cyl_cell_mass(field,data):
    return data["dens"]*data["cyl_cell_volume"]


def _cyl_radial_velocity(field,data):
    r = data["r"]
    z = data["z"]
    velx = data["velx"]
    vely = data["vely"]
    phi = np.arctan(z/r)
    velr = np.cos(phi)*velx + np.sin(phi)*vely
    return velr

def _cyl_tangential_velocity(field,data):
    r = data["r"]
    z = data["z"]
    velx = data["velx"]
    vely = data["vely"]
    phi = np.arctan(z/r)
    tanr = -np.sin(phi)*velx + np.cos(phi)*vely
    return tanr


# no need for 3D cartesian


def add_sph_fields(ds):
    """
    add yt fields for 1D spherical coordinates
    """
    ds.add_field("radial_velocity",function=_sph_radial_velocity,units='cm/s')
    ds.add_field("tangential_velocity",function=_sph_tangential_velocity,units='cm/s')
    ds.add_field("radius",function=_sph_radius,units='cm')
    ds.add_field("sph_radius",function=_sph_radius,units='cm')
    ds.add_field('sph_cell_volume',function=_sph_volume,units='cm**3')
    ds.add_field('sph_cell_mass',function=_sph_cell_mass,units='g')
    return ds

def add_cyl_fields(ds):
    """
    add yt fields for 2D cylindrical coordinates
    """
    ds.add_field("radial_velocity",function=_cyl_radial_velocity,units='cm/s')
    ds.add_field("tangential_velocity",function=_cyl_tangential_velocity,units='cm/s')
    ds.add_field("radius",function=_cyl_radius,units='cm')
    ds.add_field("cyl_radius",function=_cyl_radius,units='cm')
    ds.add_field('cyl_cell_volume',function=_cyl_volume,units='cm**3')
    ds.add_field('cyl_cell_mass',function=_cyl_cell_mass,units='g')
    return ds

def add_car_fields(ds):
    """
    add yt fields for 3D Cartesian fields coordinates
    """
    return ds

def add_ccsn_fields(ds,dim):
    if dim==1:
        ds = add_sph_fields(ds)
    elif dim==2:
        ds = add_cyl_fields(ds)
    elif dim==3:
        ds = add_car_fields(ds)
    else:
        print("Error: no such dimensions.",dim)
        quit()
    return ds

if __name__=='__main__':

    from .utility.timer import *


    # TEST 1D
    #dim   = 1
    #path  = '/Volumes/Fomalhaut-01/runs/'
    #key   = 'ccsn1d/170322_s20GR_iter/output'
    #fname = '/ccsn1d_hdf5_chk_0350'

    # TEST 2D
    #dim   = 2
    #path  = '/Volumes/Fomalhaut-01/runs/'
    #key   = 'ccsn2d/170322_s20GR_iter/output'
    #fname = '/ccsn2d_hdf5_plt_cnt_0415'

    # TEST 3D
    #dim   = 3
    #path  = "/Users/pan/Documents/runs/"
    #key   = "ccsn3d/20170224_s40GR_LS220"
    #fname = '/ccsn3d_hdf5_plt_cnt_0462' 

    fn = path+key+fname

    ds = yt.load(fn)
    #ds = add_eint(ds)

    #ds = add_sph_fields(ds)
    #ds = add_cyl_fields(ds)
    #ds = add_car_fields(ds)
    ds = add_ccsn_fields(ds,dim)
    dd = ds.all_data()

    print((list(dd.quantities.keys())))
    print((ds.field_list))
    print((ds.derived_field_list))


