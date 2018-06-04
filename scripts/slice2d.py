import yt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def slice2d(ds,var,rmax,clim=None,take_log=False):
    """
    plot a 2d slice
    """
    slice2d = yt.SlicePlot(ds,"theta",var,origin=('center','left','window'))
    slc_frb = slice2d.data_source.to_frb((2*rmax,"cm"),1024,center=(0,0,0),height=(2*rmax,"cm"))
    fig = plt.figure(1,figsize=(7,10))
    plt.cla()
    plt.clf()
    if take_log==False:
        plt.imshow(slc_frb[var].d,extent=[-rmax,rmax,-rmax,rmax],
               interpolation='nearest',
               aspect=1.0,
               cmap='Spectral_r',
               origin='center')
    else:
        plt.imshow(slc_frb[var].d,extent=[-rmax,rmax,-rmax,rmax],
               interpolation='nearest',
               aspect=1.0,
               cmap='Spectral_r',
               norm=mpl.colors.LogNorm(),
               origin='center')
            
    plt.xlim([0,rmax])
    plt.ylim([-rmax,rmax])
    if clim != None:
        plt.clim(clim)
    cbar = plt.colorbar(shrink=1.0)
    cbar.ax.set_ylabel(var,rotation=270,labelpad=25)
    plt.xlabel("R [cm]")
    plt.ylabel("Z [cm]")
    return fig
