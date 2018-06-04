#!/Users/pan/anaconda/envs/python3/bin/python
import os, sys
import yt
import numpy as np
import string
import matplotlib as mpl
import matplotlib.pyplot as plt
from optparse import OptionParser
"""
This version only supports data in 1D spherical, 2D cylindrical or 3D cartesian coordinates

Kuo-Chuan Pan  
2018.06.04
"""
def default(str):
    return str + ' [Default: %default]'
def readCommand(argv):
    usageStr = """

    -----------------------------------------------------
    A lazy script to quickly visualize FLASH data with yt

    version 0.2

    -----------------------------------------------------
    USAGE: ./yt_slice.py -n <file name> <options>

    EXAMPLE: ./yt_slice.py -n ccsn2d_hdf5_plt_cnt_0200 -v entr --rmax=3e7

    """
    parser = OptionParser(usageStr)
    parser.add_option('-n','--fname',dest="fname",
            help=default('File name'),default='')
    parser.add_option('-v','--var',dest="var",
            help=default('Plot variable.'),default='dens')
    parser.add_option('-r','--rmax',dest="rmax",
            help=default('Max radius to plot'),default=4e7)
    parser.add_option('-l','--log',dest="log",
            help=default('In log scale'),default="None")

    options, otherjunk = parser.parse_args(argv)
    if len(otherjunk) != 0:
        raise Exception('Command line input not understand: '+str(otherjunk))

    return options

def draw_slice():

    options = readCommand(sys.argv[1:])
    fn = options.fname
    ds = yt.load(fn)

    # find the dimension of the data
    dim_raw = ds.domain_dimensions
    dim_x = dim_raw[0]
    dim_y = dim_raw[1]
    dim_z = dim_raw[2]
    dim   = 0
    if dim_z == 1:
        if dim_y == 1:
            dim = 1
        else:
            dim = 2
    else:
        dim = 3
    
    if dim ==1:
        draw_1d(ds,options)
    elif dim ==2:
        draw_2d(ds,options)
    elif dim==3:
        draw_3d(ds,options)
    else:
        print("Dimension error. ", dim)
        quit()

    print("Done.")
    return

def draw_1d(ds,options):

    rmax = float(options.rmax)
    var  = options.var
    log  = options.log

    ray = ds.ray([0,0,0],[rmax,0,0])
    plt.figure()
    plt.plot(ray['t']*rmax,ray[var],'-')
    plt.xlabel("Radius [cm]")
    plt.ylabel(var)
    if log != "None":
        plt.yscale('log')
    plt.show()

    return
def draw_2d(ds,options):

    rmax = float(options.rmax)
    var  = options.var
    log  = options.log
    clim = "auto"
    
    pc = yt.SlicePlot(ds,"theta",var,origin=('center','left', 'window'))
    slc_frb = pc.data_source.to_frb((2*rmax,"cm"),1024,
            center=(0,0),height=(2*rmax,"cm"))
    plt.figure(1,figsize=(6,8))
    if var=="deps":
        my_map = "seismic"
    else:
        my_map = "Spectral_r"
    if log=="None":
        plt.imshow(slc_frb[var].d,extent=[-rmax,rmax,-rmax,rmax],
                    interpolation='nearest',
                    aspect=1.0,
                    cmap=my_map,
                    origin='center')
    else:
        plt.imshow(slc_frb[var].d,extent=[-rmax,rmax,-rmax,rmax],
                    interpolation='nearest',
                    aspect=1.0,
                    cmap=my_map,
                    norm=mpl.colors.LogNorm(),
                    origin='center')
    cbar = plt.colorbar(shrink=0.96)
    cbar.ax.set_ylabel(var,rotation=270,labelpad=15)
    if clim != "auto":
        plt.clim([cmin,cmax])
    plt.xlabel("R [cm]")
    plt.ylabel("Z [cm]")
    plt.xlim(0,rmax)
    plt.ylim(-rmax,rmax)
    plt.tight_layout()
    plt.show()

    return
def draw_3d(ds,options):

    rmax = float(options.rmax)
    var  = options.var
    log  = options.log
    clim = "auto"
    fout = "slice.png"

    ds.periodicity = (True, True, True)
    slice = yt.SlicePlot(ds,'z',[var],
            center=[0,0,0],
            width=(rmax,'cm'))

    if log!="None":
        slice.set_log(var,True)
    else:
        slice.set_log(var,False)
    if clim != "auto":
        slice.set_zlim(var,cmin,cmax)
    slice.set_cmap('entr',cmap="Spectral_r")
    slice.set_width((2.*rmax,2.*rmax))
    slice.save(fout)
    
    # open the image
    image = plt.imread(fout)
    fig, ax = plt.subplots()
    ax.imshow(image)
    ax.axis('off')
    plt.show()

    return

if __name__=='__main__':

    draw_slice()
