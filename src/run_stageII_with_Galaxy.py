"""
   Minimalistic routine for running a gravity code
"""
import sys
from amuse.lab import *
from amuse.couple import bridge
from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

#from amuse.community.symple.interface import symple	# symplectic
from amuse.community.huayno.interface import Huayno	# symplectic

from amuse.ic import make_planets_oligarch
from amuse.ext import solarsystem 

from matplotlib import pyplot as plt
import numpy as np
import math

from run_cluster import MilkyWay_galaxy
from run_solar_system_in_Galaxy import PlanetarySystemIntegrationWithPerturbers


def integrate_planetary_system_LonelyPlanets(fperturbers, Nnn, Nasteroids,
                                             time_end):

    cluster_code = PlanetarySystemIntegrationWithPerturbers(nperturbers=Nnn)
    cluster_code.read_perturber_list(fperturbers)
    cluster_code.add_planetary_system(Nasteroids)

    include_stellar_evolution = False
    if include_stellar_evolution:
        cluster_code.start_stellar_code()
    cluster_code.start_gravity_code()

    cluster_code.determine_orbital_parameters()

    galaxy_code = MilkyWay_galaxy()
    
    gravity = bridge.Bridge(verbose=False, use_threading=False)
    gravity.add_system(cluster_code, (galaxy_code,))
    gravity.timestep=0.1|units.Myr
    
    model_time = 0|units.Myr
    cluster_code.write_planetary_system()
    
    dt_diag = 0.1 |units.Myr
    t_diag = dt_diag
    dt = 0.1|units.Myr
    #time_end = cluster_code.get_last_perturber_time()
    #print("End time=", time_end.in_(units.Myr))
    while model_time<time_end:
        model_time += dt

        gravity.evolve_model(model_time)

        if model_time>t_diag:
            t_diag += dt_diag
            cluster_code.write_planetary_system()
        sys.stdout.flush()
            
    gravity.stop()


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Nnn", dest="Nnn", type="int",default = -1,
                      help="number of nearest neighbors [%default]")
    result.add_option("--Nast", dest="Nasteroids",
                      type="int",default = 10,
                      help="number of asteroids [%default]")
    result.add_option("-f", 
                      dest="fperturbers",
                      default = "lps_key_13544424148568912489.amuse",
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

    np.random.seed(31415)
    integrate_planetary_system_LonelyPlanets(o.fperturbers,
                                             o.Nnn,
                                             o.Nasteroids,
                                             o.time_end)
    
