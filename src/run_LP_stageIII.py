"""
   Minimalistic routine for running a gravity code
"""
import os
import sys
import glob
from amuse.lab import *
#from amuse.community.petar.interface import Petar
from amuse.couple import bridge
from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

#from amuse.community.symple.interface import symple	# symplectic
from amuse.community.huayno.interface import Huayno	# symplectic

from amuse.ext.orbital_elements import get_orbital_elements_from_arrays
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ic import make_planets_oligarch
from amuse.ext import solarsystem 

import numpy as np
import math

from run_LP_stageI import MilkyWay_galaxy

def get_escapers_at_time(escaper_data, model_time):
    escapers = Particles()
    for all_escapers in escaper_data:
        for pi in reversed(list(all_escapers.iter_history())):
            print("Escape time:", pi[0].age.in_(units.Myr), "N=", len(pi))
            print("at age=", pi[0].age.in_(units.Myr), model_time.in_(units.Myr))
            if pi[0].age==model_time:
                escapers.add_particles(pi)
                break
    print("escapers:", escapers)
    return escapers

        
def run_stage_III(cluster_file_names,
                  escaper_file_names,
                  time_end=10 | units.Myr,
                  dt=0.001|units.Myr):

    galaxy_code = MilkyWay_galaxy()

    print("file=", cluster_file_names)
    cluster_files = glob.glob(cluster_file_names)
    cluster_files = np.sort(cluster_files)
    print(cluster_files)

    escaper_files = glob.glob(escaper_file_names)
    print(escaper_files)
    escaper_data = []
    for escaper_file in escaper_files:
        escaper_data.append(read_set_from_file(escaper_file, close_file=True))

    converter=nbody_system.nbody_to_si(1000|units.MSun, 1|units.pc)

    cluster_file_index = 0    
    cluster_particles = read_set_from_file(cluster_files[cluster_file_index], close_file=True)

    #cluster_code = ph4(converter, number_of_workers=12)
    cluster_code = Huayno(converter)
    cluster_code.particles.add_particles(cluster_particles)

    channel_from_gd = cluster_code.particles.new_channel_to(cluster_particles)
    channel_to_gd = cluster_particles.new_channel_to(cluster_code.particles)
    
    gravity = bridge.Bridge(verbose=False, use_threading=False)
    gravity.add_system(cluster_code, (galaxy_code,))
    gravity.timestep=0.01|units.Myr

    cluster_file_index += 1
    new_cluster_particles = read_set_from_file(cluster_files[cluster_file_index],close_file=True)

    t_diag = 0|units.Myr
    dt_diag = 1 | units.Myr
    model_time = 0 | units.Myr
    dt = 0.01 |units.Myr
    while model_time < time_end:
        print("-model time=", model_time.in_(units.Myr), cluster_particles[0].age.in_(units.Myr), new_cluster_particles[0].age.in_(units.Myr))
        while model_time<new_cluster_particles[0].age:
            escaping_particles = get_escapers_at_time(escaper_data,
                                                      model_time-dt)
            #-2dt offset brings the escaper closer to the star?
            model_time += dt
            if len(escaping_particles)>0:
                s = cluster_particles[cluster_particles.key==int(17472204436860659206)]
                print("model time=", model_time.in_(units.Myr), s.age.in_(units.Myr), escaping_particles.age.in_(units.Myr), cluster_particles[0].age.in_(units.Myr), new_cluster_particles[0].age.in_(units.Myr))
                print("Host star:", s)
                print("Host star pos:", s.position.in_(units.pc))
                print("escapers pos:", escaping_particles.position.in_(units.pc))

                print("Add escapers to cluster code.")
                cluster_particles.add_particles(escaping_particles)
                cluster_particles.synchronize_to(cluster_code.particles)
            else:
                print("No escapers added.")

            E0 = cluster_code.kinetic_energy + cluster_code.potential_energy
            print("model time=", model_time.in_(units.Myr))
            gravity.evolve_model(model_time)
            channel_from_gd.copy()
            cluster_particles.age = model_time
            print("pos=", cluster_particles[0].position.in_(units.pc))
            
            N = len(cluster_code.particles)
            print(f"Evolve to time = {(model_time+dt).in_(units.Myr)}, N={N}")
            
            E1 = cluster_code.kinetic_energy + cluster_code.potential_energy
            dE = E0-E1
            ddE = (E0-E1)/E0

            print("At time=", model_time.in_(units.Myr), "dE", dE, ddE)
            sys.stdout.flush()

        print("check for new cluster model file. at t=", model_time.in_(units.Myr))
        fcluster = "freefloaters_i{0:06}.amuse".format(cluster_file_index)
        write_set_to_file(cluster_particles,
                          fcluster,
                          close_file=True)

        cluster_file_index += 1
        if len(cluster_files)<=cluster_file_index:
            print("End of cluster file reached")
            break
        cluster_file = cluster_files[cluster_file_index]
        print("pos check:", cluster_particles[0].x.in_(units.pc))
        print("pos check:", new_cluster_particles[0].x.in_(units.pc))
        new_cluster_particles.new_channel_to(cluster_particles).copy()
        new_cluster_particles = read_set_from_file(cluster_file,
                                                   close_file=True)
            
    gravity.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float",
                      default = 1|units.Myr,
                      help="end time [%default]")
    result.add_option("--dt", unit=units.Myr,
                      dest="dt", type="float",default = 0.01|units.Myr,
                      help="time step [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    np.random.seed(31415)
    run_stage_III(cluster_file_names="../S1/cluster_i*.amuse",
                  escaper_file_names="../S2/lp_escapers_key_*.amuse",
                  time_end=o.t_end,
                  dt=o.dt)
    

