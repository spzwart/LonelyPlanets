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

from amuse.ext.orbital_elements import get_orbital_elements_from_arrays
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ic import make_planets_oligarch
from amuse.ext import solarsystem 

from matplotlib import pyplot as plt
import numpy as np
import math

def semi_to_orbital_period(a, Mtot) :
    return 2*math.pi * (a**3/(constants.G*Mtot)).sqrt()

def get_orbital_elements_of_planetary_system(star, planets):
    total_masses = planets.mass + star.mass
    rel_pos = planets.position-star.position
    rel_vel = planets.velocity-star.velocity
    sma, ecc, true_anomaly,\
        inc, long_asc_node, arg_per_mat =\
            get_orbital_elements_from_arrays(rel_pos,
                                             rel_vel,
                                             total_masses,
                                             G=constants.G)
    planets.semimajor_axis = sma
    planets.eccentricity = ecc
    planets = planets[planets.eccentricity>=1]
    planets.eccentricity = 1
    planets.semimajor_axis = 1.e+10 | units.au
    planets = planets[np.isnan(planets.eccentricity)]
    planets.eccentricity = 1
    planets.semimajor_axis = 1.e+10 | units.au


def write_planetary_system(model_time, planetary_system):
    planetary_system.age = model_time
    star = planetary_system[planetary_system.name=="Sun"][0]
    filename = "iso_planets_key_"+str(star.key)+".amuse"
    #planets = planetary_system-star
    print(f"Write file {filename}") 
    write_set_to_file(planetary_system,
                      filename,
                      close_file=True,
                      append_to_file=True)
    print("File written:", filename, "N_planets=", len(planetary_system))

def remove_lost_planets(gravity, planetary_system):

    model_time = gravity.model_time
    
    sun = planetary_system[planetary_system.name=="Sun"][0]
    panda = planetary_system-sun
    asteroids = panda[panda.name=="asteroid"]
    planets = panda - asteroids
    get_orbital_elements_of_planetary_system(sun, planets)
    get_orbital_elements_of_planetary_system(sun, asteroids)

    lost_planets = panda[panda.eccentricity>1]
    la = Particles()
    if len(lost_planets)>0:
        p = lost_planets[lost_planets.type=="planet"]
        a = lost_planets[lost_planets.type=="asteroid"]
        print(f"At time={model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) unbound (planets, asteroids)")
        
        #print(lost_planets)
        lp = lost_planets.get_intersecting_subset_in(gravity.particles)
        la.add_particles(lp[lp.name=="asteroid"])
        gravity.particles.remove_particles(lp)
        planetary_system.remove_particles(lost_planets)

    lost_planets = panda[np.isnan(panda.eccentricity)]
    if len(lost_planets)>0:
        p = lost_planets[lost_planets.type=="planet"]
        a = lost_planets[lost_planets.type=="asteroid"]
        print(f"At time={model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) nan (planets, asteroids)")
        
        #print(lost_planets)
        lp = lost_planets.get_intersecting_subset_in(gravity.particles)
        la.add_particles(lp[lp.name=="asteroid"])
        gravity.particles.remove_particles(lp)
        planetary_system.remove_particles(lost_planets)


    merged_planets = Particles()
    for pi in panda:
        #print("pi=", pi.eccentricity)
        if pi.semimajor_axis*(1-pi.eccentricity)<100|units.RSun:
            merged_planets.add_particle(pi)

    if len(merged_planets)>0:
        p = merged_planets[merged_planets.type=="planet"]
        a = merged_planets[merged_planets.type=="asteroid"]
        print(f"At time={model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) merged (planets, asteroids)")
        
        mp = merged_planets.get_intersecting_subset_in(gravity.particles)
        gravity.particles.remove_particles(mp)
        planetary_system.remove_particles(merged_planets)

def integrate_planetary_system(Nasteroids):

    star = Particle()
    star.type = "star"
    star.name = "Sun"
    star.mass = 1|units.MSun
    star.position = (0,0,0) | units.au
    star.velocity = (0,0,0) | units.kms
    
    planetary_system = generate_planetary_system(star.as_set(),
                                                 Nasteroids)

    minor_bodies = planetary_system-star
    planets = minor_bodies[minor_bodies.type=="planet"]
    asteroids = minor_bodies-planets
    get_orbital_elements_of_planetary_system(star, planets)
    print(planets.semimajor_axis.in_(units.au))
    sun = planetary_system[planetary_system.name=="Sun"][0]
    
    stellar = SeBa()
    stellar.particles.add_particles(sun.as_set())
    channel_from_se = stellar.particles.new_channel_to(planetary_system)
    channel_to_se = planetary_system.new_channel_to(stellar.particles)

    stellar.evolve_model(0|units.Myr)
    channel_from_se.copy()
    
    Porb = semi_to_orbital_period(1|units.au, star.mass.sum())
    converter=nbody_system.nbody_to_si(star.mass.sum(), Porb)
    
    #gravity = symple(converter)
    #gravity = Huayno(converter, number_of_workers=8)
    gravity = Huayno(converter, mode="openmp")
    gravity.parameters.timestep_parameter = 0.03
    gravity.epsilon_squared = (0.1|units.au)**2
    print(gravity.parameters)
    gravity.particles.add_particles(planetary_system)

    #for symple
    #gravity.parameters.integrator = 10
    #gravity.parameters.timestep_parameter = 0.

    #collision_detection = gravity.stopping_conditions.collision_detection
    #collision_detection.enable()
    
    channel_from_gd = gravity.particles.new_channel_to(planetary_system)
    channel_to_gd = planetary_system.new_channel_to(gravity.particles)

    model_time = 0|units.Myr
    write_planetary_system(model_time, planetary_system)
    
    dt_diag = 1 |units.Myr
    t_diag = dt_diag
    t_end = 100|units.Myr
    dt = 0.01 | units.Myr
    while model_time<t_end:
        model_time += dt

        stellar.evolve_model(model_time)
        channel_from_se.copy()
        channel_to_gd.copy_attributes(["mass"])

        print(f"Time={model_time.in_(units.Myr)}: N= {len(gravity.particles)} = {len(planetary_system)}")

        gravity.evolve_model(model_time)
        channel_from_gd.copy()
        
        remove_lost_planets(gravity, planetary_system)
        
        print(f"Time={model_time.in_(units.Myr)}: a={planets.semimajor_axis.in_(units.au)}, e={planets.eccentricity}")
        if model_time>t_diag:
            t_diag += dt_diag
            write_planetary_system(model_time, planetary_system)
        sys.stdout.flush()
            
    gravity.stop()

def generate_planetary_system(parent_star, Nasteroids):

    planetary_system = new_solar_system()
    sun = planetary_system[planetary_system.name=="SUN"][0]
    parent_star.mass = sun.mass
    planetary_system.position-=sun.position
    planetary_system.velocity-=sun.velocity
    planetary_system.remove_particle(sun)
    #planets = planetary_system[2:]
    planets = planetary_system[4:]
    planets = planets[:-1]
    planets.type = "planet"
    planets.position += parent_star.position
    planets.velocity += parent_star.velocity

    converter=nbody_system.nbody_to_si(parent_star.mass.sum(), 1|units.au)
    asteroids = ProtoPlanetaryDisk(Nasteroids,
                                   densitypower=1.5, Rmin=1, Rmax=100,
                                   q_out=1, discfraction=0.01,
                                   convert_nbody=converter).result
    asteroids.mass = 0 | units.MEarth
    asteroids.name = "asteroid"
    asteroids.type = "asteroid"
    asteroids.position += parent_star.position
    asteroids.velocity += parent_star.velocity
    
    parent_star.add_particles(planets)
    parent_star.add_particles(asteroids)
    return parent_star

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Nasteroids", type="int",default = 1000,
                      help="number of asteroids [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    np.random.seed(31415)
    integrate_planetary_system(o.Nasteroids)
    
