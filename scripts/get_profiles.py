import yt
import numpy as np
from . import add_fields as af
from scipy.interpolate import interp1d


class RadialProfile():
    """
    Get uniform spaced radial profiles 
    Note: here we assume center located at [0,0,0]
    """
    def __init__(self,dr=1.e5,rmax=5.e7):
        self.dr     = dr
        self.rmax   = rmax
        self.nbins  = int(rmax/dr)
        self.radius = np.linspace(dr,rmax,self.nbins)
        self.profiles = {}
        return

    def get_1d_profile(self,ds,variables):
        """
        get 1d profile

        input: variables: [str]
        """
        rmax = self.rmax
        profiles = self.profiles
        ctr  = [0,0,0]
        pt   = [rmax,0,0]
        ray = ds.ray(ctr,pt)
        raw_radius = (ray['t'].in_cgs().v)*rmax
        s1 = raw_radius.argsort()
        raw = {}
        for var in variables:
            raw[var] = ray[var].in_cgs().v
        radius = np.array(raw_radius)[s1]
        for var in variables:
            value  = np.array(raw[var][s1])
            f1d = interp1d(radius,value)
            profiles[var] = np.zeros(len(self.radius))
            for i,r in enumerate(self.radius):
                try:
                    profiles[var][i] = f1d(r)
                except:
                    if r <= radius[0]:
                        profiles[var][i]=value[0]
                    if r > radius[-1]:
                        profiles[var][i]=value[-1]

        self.profiles = profiles
        return

    def get_2d_profile(self,ds,variables):
        """
        get radial profiles from 2d data
        """
        source = ds.sphere([0,0,0], (self.rmax,"cm"))
        yt_profile = yt.create_profile(source,
                "cyl_radius",
                variables,
                n_bins=self.nbins,
                extrema={'cyl_radius':(self.dr/2,self.rmax+self.dr/2)},
                logs={'cyl_radius':False},
                weight_field='cyl_cell_mass',
                accumulation=False)
        for var in variables:
            self.profiles[var]=yt_profile[var].v
        return

    def get_3d_profile(self,ds,variables):
        """
        get radial profiles from 3d data
        """
        source = ds.sphere([0,0,0], (self.rmax,"cm"))
        yt_profile = yt.create_profile(source,
                "radius",
                variables,
                n_bins=self.nbins,
                extrema={'radius':(self.dr/2,self.rmax+self.dr/2)},
                logs={'radius':False},
                weight_field='cell_mass',
                accumulation=False)
        for var in variables:
            self.profiles[var]=yt_profile[var].v
        return

    def get_profile(self,ds,dim,variables):
        if dim==1:
            self.get_1d_profile(ds,variables)
        elif dim==2:
            self.get_2d_profile(ds,variables)
        elif dim==3:
            self.get_3d_profile(ds,variables)
        else:
            print("Error: no such dimension.", dim)
            quit()
        return

if __name__=='__main__':

    import matplotlib.pyplot as plt
    from .eos import NuclearEos as ne
    eos_file = "/Volumes/Fomalhaut-01/eos/SFHo.h5"
    neos     = ne.NuclearEOS(eos_file)

    variables = ["dens","entr", "ye  ","deps"]

    test1D = False
    if test1D:
        print(" *** unit test: 1D ***")
        path = "/Volumes/Fomalhaut-01/runs/ccsn1d/170322_s20GR_iter/output"
        fn   = path+"/ccsn1d_hdf5_chk_0070"
        ds = yt.load(fn)
        af.neos = neos
        ds = af.add_ccsn_fields(ds,1)

        rp = RadialProfile()
        rp.get_profile(ds,1,variables)
        plt.figure(1)
        plt.plot(rp.radius,rp.profiles["dens"])
        plt.plot(rp.radius,rp.profiles["entr"])
        plt.plot(rp.radius,rp.profiles["ye  "])
        plt.yscale('log')
        plt.xlim([0,2e7])
        plt.show()

    test2D = False
    if test2D:
        print(" *** unit test: 2D ***")
        path = "/Volumes/Fomalhaut-01/runs/ccsn2d/170322_s20GR_iter/output"
        fn   = path+"/ccsn2d_hdf5_plt_cnt_0070"
        ds = yt.load(fn)
        af.neos = neos
        ds = af.add_eint(ds)
        ds = af.add_ccsn_fields(ds,2)
    
        rp = RadialProfile()
        rp.get_profile(ds,2,variables)
        plt.figure(2)
        plt.plot(rp.radius,rp.profiles["dens"])
        plt.plot(rp.radius,rp.profiles["entr"])
        plt.plot(rp.radius,rp.profiles["ye  "])
        plt.yscale('log')
        plt.xlim([0,2e7])
        plt.show()

    test3D = True
    if test3D:
        print(" *** unit test: 3D ***")
        path = "/Users/pan/Documents/runs/ccsn3d/20170224_s40GR_LS220"
        fn   = path+"/ccsn3d_hdf5_plt_cnt_0462"
        ds = yt.load(fn)
        af.neos = neos
        ds = af.add_ccsn_fields(ds,3)
    
        rp = RadialProfile()
        rp.get_profile(ds,3,variables)
        plt.figure(3)
        #plt.plot(rp.radius,rp.profiles["dens"])
        plt.plot(rp.radius,rp.profiles["entr"])
        plt.plot(rp.radius,rp.profiles["ye  "])
        #plt.yscale('log')
        plt.xlim([0,5e7])
        plt.show()
