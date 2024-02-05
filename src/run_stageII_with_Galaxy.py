"""
   Minimalistic routine for running a gravity code

Huayno integrator types
CONSTANT = 0, SHARED2 = 1, PASS_KDK = 2, HOLD_KDK = 3, BRIDGE_KDK = 4, 
      EXTRAPOLATE = 5, PASS_DKD = 7, HOLD_DKD = 8, PPASS_DKD = 9, BRIDGE_DKD = 10,
      CC = 11, CC_KEPLER = 12, OK = 13, KEPLER = 14, SHARED4 = 15, FOURTH_M4 = 16, FOURTH_M5 = 17,
      SHARED6 = 18, SHARED8 = 19, SHARED10 = 20, SHAREDBS = 21, CCC = 22, CCC_KEPLER = 23,
      CC_BS = 24, CCC_BS = 25, BS_CC_KEPLER = 26, CC_BSA = 27, CCC_BSA = 28, SHARED2_COLLISIONS = 29,
      SHARED4_COLLISIONS = 30, SHARED6_COLLISIONS = 31, SHARED8_COLLISIONS = 32, 
      SHARED10_COLLISIONS = 33, CONSTANT2 = 34, CONSTANT4 = 35, CONSTANT6 = 36, 
      CONSTANT8 = 37, CONSTANT10 = 38, ERROR_CONTROL=39, CC_SHARED10=40, CCC_SHARED10=41
"""
import sys
import os
import os.path
import time as wallclock
from amuse.lab import *
from amuse.couple import bridge
from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

#from amuse.community.symple.interface import symple	# symplectic
from amuse.community.huayno.interface import Huayno	# symplectic

from amuse.ic import make_planets_oligarch
from amuse.ext import solarsystem 

#from matplotlib import pyplot as plt
import numpy as np
import math

from run_cluster import MilkyWay_galaxy
from run_solar_system_in_Galaxy import PlanetarySystemIntegrationWithPerturbers


def integrate_planetary_system_LonelyPlanets(fperturbers, Nnn, Nasteroids,
                                             time_end, integrator):

    wct_initialization = wallclock.time()    
    dt_diag = 1.0 | units.Myr
    cluster_code = PlanetarySystemIntegrationWithPerturbers(nperturbers=Nnn,
                                                            maximal_timestep=dt_diag)
    cluster_code.read_perturber_list(fperturbers)
    cluster_code.add_planetary_system(Nasteroids)

    include_stellar_evolution = False
    if include_stellar_evolution:
        cluster_code.start_stellar_code()
    cluster_code.start_gravity_code(integrator)

    cluster_code.determine_orbital_parameters()

    #galaxy_code = MilkyWay_galaxy(Mb=0*1.40592e10| units.MSun,
    #                              Md=0*8.5608e10| units.MSun,
    #                              Mh=0*1.07068e11 | units.MSun)
    galaxy_code = MilkyWay_galaxy()
    
    gravity = bridge.Bridge(verbose=False, use_threading=False)
    gravity.add_system(cluster_code, (galaxy_code,))
    gravity.timestep=0.5*dt_diag

    cluster_code.wct_initialization = (wallclock.time()-wct_initialization) | units.s    
    model_time = 0|units.Myr
    cluster_code.write_planetary_system()
    
    dt = min(time_end, dt_diag)
    t_diag = dt_diag
    while model_time<time_end:
        model_time += dt

        gravity.evolve_model(model_time)

        if model_time>t_diag:
            t_diag += dt_diag
            cluster_code.write_planetary_system()
        sys.stdout.flush()

        if cluster_code.terminating_criterium_reached():
            print("Stop the run due to lost planets/asteroids.")
            time_end = t_diag


        if os.path.isfile("STOP"):
            file = open("STOP", "r")
            keys = file.readlines()
            file.close() 
            for key in keys:
                if str(key) in str(cluster_code.key):
                    print("Stop the run by manual intervention.")
                    time_end = t_diag
            #os.remove("STOP") 

    gravity.stop()
    cluster_code.print_wallclock_time()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Nnn", dest="Nnn", type="int",default = -1,
                      help="number of nearest neighbors [%default]")
    result.add_option("--Nast", dest="Nasteroids",
                      type="int",default = 1000,
                      help="number of asteroids [%default]")
    result.add_option("--seed", dest="seed",
                      type="int",default = 12345,
                      help="random number seed [%default]")
    result.add_option("-I", dest="integrator",
                      type="int", default = 8,
                      help="integrator id (Huayno) [%default]")
    result.add_option("-f", 
                      dest="fperturbers",
                      #default = "lps_key_13544424148568912489.amuse",
                      default = "",
                      help="stellar input filename [%default]")
    result.add_option("-t", 
                      dest="time_end",
                      unit=units.Myr,
                      type="float",
                      default = 1|units.Myr,
                      help="Simulation end time [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    np.random.seed(o.seed)
    integrate_planetary_system_LonelyPlanets(o.fperturbers,
                                             o.Nnn,
                                             o.Nasteroids,
                                             o.time_end,
                                             o.integrator)
    
