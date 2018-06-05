import matplotlib as mpl
mpl.use('Agg')
import yt
yt.enable_parallelism()
import string
import numpy as np
from yt.visualization.volume_rendering.transfer_function_helper import *
from yt.visualization.volume_rendering.api import LineSource
from yt.utilities.amr_kdtree.api import AMRKDTree
from yt.units.dimensions import mass, energy, temperature
from yt.units import cm
import sys
sys.path.insert(0,'..')
from my_volume_rendering_setting import *
"""
Volume rendering plot.

Change sim_volume_rendering_setting.py for your simulation

"""

def get_fn(path,header,cycle):
    """
    return the file name. here we assume using plt files
    """
    return path+'/'+header+'_hdf5_plt_cnt_'+string.zfill(cycle,4)

def plot_a_vr(path,header,cycle,use_ghost_zones=True,annotate=True,rotate=False,zoom=False):
    """
    volume rendering plot

    path [string]  : the path of data files
    header [string]: the header of the data file
    cycle [int]    : the data cycle

    ex. for file: /data/ccsn3d_hdf5_plt_cnt_0100

        path   = "/data"
        header = "ccsn3d"
        cycle  = 100
    """

    fn = get_fn(path,header,cycle)
    ds = yt.load(fn)

    # register new units for entropy
    ds.unit_registry.add('kB',1.3806488e-16,dimensions=energy/temperature,tex_repr='k_{B}')
    ds.unit_registry.add('by',1.674e-24,dimensions=mass,tex_repr='baryon')

    # create a new derived field: Entropy
    def _entr(field,data):
        """
	fixed the entropy units in flash
        """
        entr = data["entr"]
        kb_by = yt.YTQuantity(1,'erg')/yt.YTQuantity(1,'K')/yt.YTQuantity(1,'g')*(1.3806488e-16/1.674e-24)
        return entr*kb_by
    ds.add_field("Entr",function=_entr,units="kB/by",
            display_name="Entropy",
            dimensions=energy/temperature/mass)

    # create a new derived field: Entropy
    def _entrdens(field,data):
        """
        create a new derived field to show both entropy and density

        if density > PNS_DENSITY:
            entropy = PNS_ENTR
        else:
            entropy = entropy

        """
        dens = data["dens"]
        entr = data["entr"]
        entrdens = entr*(np.exp(-(dens.in_cgs()/PNS_DENSITY)**5))+PNS_ENTR
        kb_by = yt.YTQuantity(1,'erg')/yt.YTQuantity(1,'K')/yt.YTQuantity(1,'g')*(1.3806488e-16/1.674e-24)
        return entrdens*kb_by
    ds.add_field("Entropy",function=_entrdens,units="kB/by",
            display_name="Entropy",
            dimensions=energy/temperature/mass)

    #debug: also plot a entropy slice 
    if yt.is_root():
        pc = yt.SlicePlot(ds,'z','Entr')
        pc.zoom(40)
        pc.set_log('Entr',False)
        pc.save('fig_slice_z_'+header+'_'+string.zfill(cycle,4)+'.png')

    # get the entropy range from time
    time = ds.current_time.in_cgs().v - TSHIFT

    if CUSTOM_EMAX:
        emin, emax = get_emin_emax(time)
    else:
        entropy_max, pe = ds.find_max('Entropy')
        emax = entropy_max.v - 1.0
        emin = emax - 3.0

    if yt.is_root():
        print "emin/emax:",emin,emax
        # check if emax > emin
        if emax <= emin:
            print "Error: emax <= emin"
            quit()

    # only render the region with r < 1.e8 cm
    # this is necessary if we want to include ghost zones
    sphere = ds.sphere([0,0,0],(1.e8, 'cm'))
    sc = yt.create_scene(sphere, field='Entropy')

    # set the camera resolution and width
    sc.camera.north_vector = [0,0,1]
    sc.camera.resolution = (VR_RESOLUTION,VR_RESOLUTION)
    sc.camera.set_width(ds.quan(VR_WIDTH,'cm'))
    #sc.camera.set_width(ds.quan(2.5e7,'cm'))

    # define the tranfer function for high entropy region
    def linearFunc(values,minv,maxv):
        try:
            return ((values) - values.min())/(values.max()-values.min())
            #return (na.sqrt(values) - values.min())/(values.max()-values.min())
        except:
            return 1.0

    source = sc.get_source(0)
    source.set_use_ghost_zones(use_ghost_zones) # *** SLOW ***
    source.set_log(False)
    source.grey_opacity=True

    # create a transfer function helper
    tfh = TransferFunctionHelper(ds)
    tfh.set_field('Entropy')
    tfh.set_bounds()
    tfh.set_log(False)
    tfh.build_transfer_function()

    # add a thin layer to show the shock front
    esh, ew = get_shock_entropy(time,emin,emax)
    tfh.tf.add_layers(1,w=(ew),mi=(esh),ma=(esh+0.5),col_bounds=[4.0,(esh+1.0)],
                    alpha=10.0*esh*np.ones(1),colormap='cool_r')

    # plot the PNS at entr = PNS_ENTR
    tfh.tf.add_layers(1,w=0.5,mi=0.49,ma=0.55,col_bounds=[0.05,0.55],
                    alpha=100.0*np.ones(1),colormap='Purples')

    # map the high entropy region: version 3
    tfh.tf.map_to_colormap(
            emin,emax,
            scale=100.0,
            scale_func=linearFunc,
            colormap='autumn')

    # version 1: use many layers
    #tfh.tf.add_layers(12,w=0.05,mi=emin,ma=emax,col_bounds=[emin,emax],
    #        alpha=70*np.linspace(0.0,1.5,10),colormap='hot')

    tfh.tf.grey_opacity = True
    source.set_transfer_function(tfh.tf)
    source.grey_opacity=True

    # plot the transfer function
    #source.tfh.plot('fig_transfer_function_entr.png', profile_field='cell_mass')

    # plot volume rendering plot without annotation 
    if not annotate:
        sc.save('fig_vr_'+header+'_'+string.zfill(cycle,4)+'.png', sigma_clip=4)

    else:
        # with annotation
        sc.annotate_axes(alpha=0.8)
        #sc.annotate_scale()  # 

        # line annotation
        #annotate_width = 1.0e7  # [cm]
        #annotate_axis  = 1.5e7  # [cm]
        #annotate_at    = (annotate_axis/2.0 - annotate_width/(2.0*np.sqrt(2.0)))

        #colors   = np.array([[0.45,0.5,0.5,0.7]]) # white
        #vertices = np.array([[[annotate_at,0.,0.],[0.,annotate_at,0.]]])*annotate_width*cm
        #lines    = LineSource(vertices,colors)
        #sc.add_source(lines)

        #sc.annotate_domain(ds,color=[1,1,1,0.01])
        #text_string= "Time = %.1f (ms)" % (float(ds.current_time.to('s')*1.e3))
        text_string= "Time = %.1f (ms)" % (float(time*1.e3))
        sc.save_annotated("fig_vr_"+header+"_annotated_"+string.zfill(cycle,4)+'.png',
                sigma_clip=4.0,
                label_fmt="%.1f",
                text_annotate=[[(0.05,0.95), 
                text_string,dict(color="w", fontsize="20", horizontalalignment="left")]])

    # rotate the camera
    if DO_ROT:
        rotate = True

    if rotate:
        fstart = 1   # modify here for restart
        frames = ROT_FRAMES # total number of frames for rotation
        cam = sc.camera
        i=0
        for _ in cam.iter_rotate(2.0*np.pi, frames):
            i+=1
            if i>= fstart:
                sc.render()
                if not annotate:
                    sc.save("fig_vr_"+header+"_"+string.zfill(cycle,4)+'_rot_'+string.zfill(i,4)+'.png',
                            sigma_clip=2)
                else:
                    sc.save_annotated("fig_vr_"+header+"_annotated_"+
                        string.zfill(cycle,4)+'_rot_'+string.zfill(i,4)+'.png',
                        sigma_clip=4.0,
                        text_annotate=[[(0.05,0.95), 
                        text_string,dict(color="w", fontsize="20", horizontalalignment="left")]])

    # TODO: zoom in or zoom out
    if DO_ZOOM:
        zoom = True

    if zoom:
        fstart = 0 # modify here for restart
        frames = ZOOM_FRAMES
        cam = sc.camera
        i=0
        for _ in cam.iter_zoom(ZOOM_FACT,ZOOM_FRAMES):
            i+=1
            if i>= fstart:
                sc.render()
                if not annotate:
                    sc.save("fig_vr_"+header+"_"+string.zfill(cycle,4)+'_zoom_'+string.zfill(i,4)+'.png',
                            sigma_clip=2)
                else:
                    sc.save_annotated("fig_vr_"+header+"_annotated_"+
                        string.zfill(cycle,4)+'_zoom_'+string.zfill(i,4)+'.png',
                        sigma_clip=4.0,
                        text_annotate=[[(0.05,0.95), 
                        text_string,dict(color="w", fontsize="20", horizontalalignment="left")]])
        quit()


    return


if __name__=='__main__':


    path   = "/Volumes/fs0/pan/runs/data"
    header = "mesa20_LR"
    cycle  = 283

    plot_a_vr(path,header,cycle,use_ghost_zones=False,annotate=False,rotate=False)
