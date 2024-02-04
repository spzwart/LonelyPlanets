"""
   Minimalistic routine for running a gravity code
"""
from amuse.lab import *
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

#from matplotlib import pyplot as plt
import numpy as np
import math

def ZAMS_radius(mass):
    log_mass = np.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun))**2
    alpha = 0.08353 + 0.0565*log_mass
    beta  = 0.01291 + 0.2226*log_mass
    gamma = 0.1151 + 0.06267*log_mass
    r_zams = pow(mass.value_in(units.MSun), 1.25) * (0.1148 + 0.8604*mass_sq) / (0.04651 + mass_sq)
    return r_zams | units.RSun

def generate_initial_star_cluster(N=10, Rvir = 1|units.pc, W0=7.0):

    Mmin = 0.08|units.MSun
    Mmax = 30|units.MSun
    mass = new_kroupa_mass_distribution(N, mass_min=Mmin, mass_max=Mmax)
    converter=nbody_system.nbody_to_si(mass.sum(), Rvir)
    
    bodies = new_king_model(N, W0, convert_nbody=converter)
    bodies.name = "star"
    bodies.type = "star"
    bodies.mass = mass
    bodies.radius = ZAMS_radius(bodies.mass)
    bodies.scale_to_standard(converter)
    bodies.move_to_center()

    return bodies

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [%default]")
    result.add_option("-W", dest="Wo", type="float",default = 5,
                      help="King model profile W0[%default]")
    result.add_option("-R", dest="Rvir", unit=units.pc,
                      type="float",default = 1|units.pc,
                      help="cluster virial radius [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    np.random.seed(31415)
    initial_bodies = generate_initial_star_cluster(o.N, o.Rvir, o.Wo)
    bodies = initial_bodies.copy()
    write_set_to_file(bodies, "initial_cluster.amuse")
