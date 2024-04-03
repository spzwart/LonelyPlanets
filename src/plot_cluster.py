"""
   Minimalistic routine for running a gravity code
"""
from amuse.lab import *
from amuse.couple import bridge
from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

from amuse.community.symple.interface import symple	# symplectic
from amuse.community.huayno.interface import Huayno	# symplectic

from amuse.ext.orbital_elements import get_orbital_elements_from_arrays
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ic import make_planets_oligarch
from amuse.ext import solarsystem 

from matplotlib import pyplot as plt
import numpy as np
import math

def plot_clock_time(fig, time_fraction):
    left, bottom, width, height = [0.8, 0.75, 0.1, 0.1]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    ax2.axis('off')

    angle = np.linspace(0, 2*np.pi, 150) 
    radius = 10|units.kpc
    x = radius * np.cos(angle) 
    y = radius * np.sin(angle) 

    time_angle = 2*np.pi*time_fraction    
    xarrow = radius * np.cos(time_angle) 
    yarrow = radius * np.sin(time_angle) 
    
    ax2.plot(x.value_in(units.kpc), y.value_in(units.kpc), lw=1, c='k')
    ax2.plot([0, xarrow.value_in(units.kpc)],
             [0, yarrow.value_in(units.kpc)], lw=1, c='k') 
    ax2.set_aspect(1) 


def plot_clock_com(fig, com):
    left, bottom, width, height = [0.8, 0.75, 0.1, 0.1]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    ax2.axis('off')

    angle = np.linspace(0, 2*np.pi, 150) 
    radius = 8.5|units.kpc
    x = radius * np.cos(angle) 
    y = radius * np.sin(angle) 
 
    ax2.plot(x.value_in(units.kpc), y.value_in(units.kpc), lw=1, c='k')
    ax2.plot([0, com.x.value_in(units.kpc)],
             [0, com.y.value_in(units.kpc)], lw=1, c='k') 
    ax2.set_aspect(1) 

def plot_perifery(fig, bodies, orbit):
    left, bottom, width, height = [0.8, 0.75, 0.15, 0.15]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    #ax2.axis('off')
    ax2.set_aspect(1) 

    c=np.log10(bodies.mass.value_in(units.MSun))
               
    s=10*(1+np.clip(np.log10(bodies.mass.value_in(units.MSun)), -1, 2))
    ec = ax2.scatter(bodies.x.value_in(units.pc),
                     bodies.y.value_in(units.pc), c=c, s=s,
                     cmap="rainbow_r",
                     vmin=-1, vmax=2)
    
    sun = bodies[bodies.name=="Sun"]
    ax2.scatter(sun.x.value_in(units.pc),
                sun.y.value_in(units.pc),
                marker=".", s=10, c='k')
    lim = 100
    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)

    #com = bodies.center_of_mass()
    #com = sun[0].position
    #com = bodies.density_center()
    converter=nbody_system.nbody_to_si(bodies.mass.sum(), 1|units.pc)
    com = densitycentre_coreradius_coredens(particles,
                                            unit_converter=converter)
    print(com)
    xxx
    
    #    bound = bound_subset(particles, tidal_radius=None, unit_converter=None, density_weighting_power=2,
    #        smoothing_length_squared=zero, G=constants.G, core=None,
    #        reuse_hop=False, hop=HopContainer(), gravity_code=None):

    print(com.in_(units.pc))
    ax2.plot([0, com.x.value_in(units.pc)],
             [0, com.y.value_in(units.pc)], lw=1, c='k', ls=":") 

    orbit = np.array(orbit).T
    x = orbit[0]
    y = orbit[1]
    ax2.plot(x, y, c='k', lw=1)

def plot_particles(bodies, alpha=1, ax=plt):
    c=np.log10(bodies.mass.value_in(units.MSun))
    s=1+10*(1+np.clip(np.log10(bodies.mass.value_in(units.MSun)), -1, 2))
    ec = ax.scatter(bodies.x.value_in(units.pc),
                    bodies.y.value_in(units.pc), c=c, s=s,
                    cmap="rainbow_r", alpha=alpha,
                    vmin=-1, vmax=2)
    return ec

def plot_cluster(bodies, orbit, index, savefig=True):

    fig, ax1 = plt.subplots(1) 
        
    #com = bodies.center_of_mass()
    sun = bodies[bodies.name=="Sun"]
    #com = sun[0].position

    from amuse.datamodel.particle_attributes import densitycentre_coreradius_coredens
    converter=nbody_system.nbody_to_si(bodies.mass.sum(), 1|units.pc)
    result = densitycentre_coreradius_coredens(bodies,
                                               unit_converter=converter)
    cod = result[0]
    rcore = result[1]
    ncore = result[2]
    print("cod=", cod.value_in(units.pc))
    print("rcore=", rcore.value_in(units.pc),
          "ncore=", ncore.in_(units.MSun/units.pc**3))
    com = cod

    rdist = (sun.position-cod).lengths()
    for ri in range(len(rdist)):
        print("stellar distance:", rdist[ri].in_(units.pc), sun[ri].key)


    co = bodies[bodies.stellar_type>=14|units.stellar_type]
    rdist = (co.position-cod).lengths()
    for ci in range(len(co)):
        print("co distance:", rdist[ci].in_(units.pc), co[ci].mass.in_(units.MSun), co[ci].key)
    
        
    check_binaries = False
    if check_binaries:
        binaries = bodies.get_binaries()
        print("N-binaries=", len(binaries))
        print(binaries)

    compact_objects = bodies[bodies.stellar_type>=14|units.stellar_type]
    stars = bodies - compact_objects
    
    from amuse.datamodel.particle_attributes import bound_subset
    bound = bound_subset(stars, tidal_radius=None,
                         unit_converter=converter,
                         density_weighting_power=2,
                         smoothing_length_squared=zero,
                         G=constants.G, core=None,
                         reuse_hop=False, 
                         gravity_code=None)
    print("bound=", len(bound))
    print("com=", bound.center_of_mass().in_(units.pc))
    unbound = stars-bound
    
    plot_clock_com(fig, com)
    #time_fraction = bodies[0].age/(220|units.Myr)
    #plot_clock(fig, time_fraction)
    if len(orbit)>0:
        plot_perifery(fig, bodies, orbit)    
    
    #bodies.move_to_center()
    #bodies.position -= sun[0].position
    #bodies.velocity -= sun[0].velocity
    #c=np.log10(bodies.mass.value_in(units.MSun))

    bodies.position -= com
    
    ec = plot_particles(bound, alpha=1, ax=ax1)
    plot_particles(unbound, alpha=1, ax=ax1)
    s=1+10*(1+np.clip(np.log10(compact_objects.mass.value_in(units.MSun)), -1, 2))
    ax1.scatter(compact_objects.x.value_in(units.pc),
               compact_objects.y.value_in(units.pc), c='k', s=s)
    

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax1)
    cax = fig.add_axes([0.82, 0.12, 0.05, 0.6])
    #plt.colorbar(im, cax=cax)
    cbar = plt.colorbar(ec, cax=cax)
    cbar.set_label('$\log_{10}(M/M_{\odot})$')

    ax1.scatter(sun.x.value_in(units.pc),
                sun.y.value_in(units.pc), c='yellow',
                marker="*", s=100, lw=1,
                edgecolors='k')

    #print(bodies[0].age.in_(units.Myr))
    print(bodies[0].position.in_(units.pc))
    escapers = bodies[bodies.type=="asteroid"]
    print(escapers.position.in_(units.pc))
    if len(escapers)>0:
        print(escapers[0].age.in_(units.Myr))
    print(escapers)
    ax1.scatter(escapers.x.value_in(units.pc),
                escapers.y.value_in(units.pc), 
                marker="o", s=10, lw=2,
                edgecolors='k', c='k')
    
    #ax1.scatter(sun[0].x.value_in(units.pc),
    #            sun[0].y.value_in(units.pc), c='yellow',
    #            marker="s", s=100, lw=1,
    #            edgecolors='k', facecolor='k')
    lim = 100
    ax1.set_xlim(-lim, lim)
    ax1.set_ylim(-lim, lim)
    ax1.set_aspect('equal')

    if savefig:
        figname = "cluster_i{0:04}.jpg".format(index)
        plt.savefig(figname)
        plt.close()
    else:
        plt.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="fcluster", default = "",
                      help="input filename [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    if len(o.fcluster)>0:
        file = o.fcluster
        index = int(file.split("_i")[1].split(".amuse")[0])
        bodies = read_set_from_file(file, close_file=True)
        orbit = []
        plot_cluster(bodies.copy(), orbit, index, savefig=False)
    else:
        from glob import glob
        files = glob("cluster_i*.amuse")
        files = np.sort(files)
        orbit = []
        for file in files:
            print(f"process {file}")
            index = int(file.split("cluster_i")[1].split(".amuse")[0])
            bodies = read_set_from_file(file, close_file=True)
            orbit.append(bodies.center_of_mass().value_in(units.pc))
            plot_cluster(bodies.copy(), orbit, index)
