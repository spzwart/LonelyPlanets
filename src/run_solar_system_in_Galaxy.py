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

from run_cluster import MilkyWay_galaxy

class PlanetarySystemIntegrationWithPerturbers(object):
    def __init__(self, nperturbers = -1):
        self.model_time = 0 | units.yr
        self.gravity_code = None
        self.stellar_code = None
        self.particles = Particles(0)
        self.nperturbers = nperturbers
        self.perturbers = Particles(0)
        self.perturber_list = None
        self.perturber_list_index = 0

        self.converter=nbody_system.nbody_to_si(1|units.MSun,
                                                1000|units.au)

    def get_perturbed_particle(self):
        perturbed_particle = self.particles[self.particles.type=="star"]
        if len(perturbed_particle)==1:
            return perturbed_particle
        else:
            print("No or too many perturbed particle(s)")
            print("Stop")
            exit(-1)

    def add_planetary_system(self, Nasteroids):
        parent_star = self.get_perturbed_particle()
        planetary_system = new_solar_system()
        sun = planetary_system[planetary_system.name=="SUN"][0]
        parent_star.mass = sun.mass
        planetary_system.position-=sun.position
        planetary_system.velocity-=sun.velocity
        planetary_system.remove_particle(sun)
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
    
        self.particles.add_particles(planets)
        self.particles.add_particles(asteroids)

    def determine_orbital_parameters(self):
        parent_star = self.particles[self.particles.type=="star"]
        minor_bodies = self.particles-parent_star
        planets = minor_bodies[minor_bodies.type=="planet"]
        asteroids = minor_bodies-planets
        get_orbital_elements_of_planetary_system(parent_star, planets)
        #print(planets.semimajor_axis.in_(units.au))
        
    def print_status(self):
        stars = self.particles[self.particles.type=="star"]
        planets = self.particles[self.particles.type=="planet"]
        asteroids = self.particles[self.particles.type=="asteroid"]
        print(f"Time={self.model_time.in_(units.Myr)}: N= {len(self.particles)} = {len(stars)}+{len(planets)}+{len(asteroids)}, Nperturbers={len(self.perturbers)}")
        
    def start_stellar_code(self, include_stellar_evolution):
        if include_stellar_evolution:
            self.stellar_code = SeBa()#redirection="none")
            stars = self.particles[self.particles.type=="star"]
            self.stellar_code.particles.add_particles(stars)
            self.from_stellar = self.stellar_code.particles.new_channel_to(self.particles)
            self.to_stellar = self.particles.new_channel_to(self.stellar_code.particles)
        else:
            self.stellar_code = None
    def start_gravity_code(self):
        self.gravity_code = Huayno(convert_nbody=self.converter,
                                   mode="openmp",
                                   redirection="none")
        self.gravity_code.particles.add_particles(self.particles)
        self.gravity_code.parameters.timestep_parameter = 0.03

        self.from_gravity = self.gravity_code.particles.new_channel_to(self.particles)
        self.to_gravity = self.particles.new_channel_to(self.gravity_code.particles)

    def get_last_perturber_time(self):
        return self.perturber_list[0][0].age
    def get_next_perturber_time(self):
        return self.perturber_list[self.perturber_list_index-1][0].age
    def get_current_perturber_time(self):
        return self.perturber_list[self.perturber_list_index][0].age

    def evolve_model(self, model_time):
        print("Evolve to", model_time.in_(units.Myr))
        
        self.model_time = self.perturber_list[self.perturber_list_index][0].age
        #print("Current model time=", self.model_time.in_(units.Myr))
        #print("Evolve to time=", model_time.in_(units.Myr))
        while self.model_time<model_time:

            print("time=", self.model_time.in_(units.Myr))

            if self.model_time<self.get_last_perturber_time():
                self.add_perturbers()

            self.print_status()

            self.evolve_system_for_one_step()

            self.remove_perturbers()
            
            ##cluster_code.remove_lost_planets()

    def evolve_system_for_one_step(self):
        print("evolve from", self.model_time.in_(units.Myr),
              "to ", self.get_next_perturber_time().in_(units.Myr),
              "index=", self.perturber_list_index)
        #print("pos=", self.model_time.in_(units.Myr), self.particles.x.in_(units.pc))
        if self.stellar_code != None:
            self.stellar_code.evolve_model(self.model_time)
            self.from_stellar.copy()
        self.to_gravity.copy()        
        self.gravity_code.evolve_model(self.get_next_perturber_time())
        self.from_gravity.copy()
        if len(self.perturbers)>0:
            self.to_perturbers.copy()
        if self.perturber_list_index==0:
            print("End of perturber list reached at time=",
                  self.model_time.in_(units.Myr))
        else:
            self.perturber_list_index -= 1

        self.model_time = self.get_current_perturber_time()
        print("Evolution done, current index=", self.perturber_list_index)
        #print("pos=", self.model_time.in_(units.Myr), self.particles.x.in_(units.pc))

    def get_parent_star(self):
        star = self.particles[self.particles.name=="Sun"]
        return star

    """
    def add_particles(self, particles):
        self.particles.add_particles(particles)
        self.particles.synchronize_to(self.gravity_code.particles)
        #self.particles.synchronize_to(self.stellar_code.particles)
    """
    def remove_particles(self, particles):
        self.particles.remove_remove(particles)
        self.particles.synchronize_to(self.gravity_code.particles)
        #self.particles.synchronize_to(self.stellar_code.particles)
    def get_gravity_at_point(self):
        return self.gravity_code.get_gravity_at_point()
    @property
    def potential_energy(self):
        return self.gravity_code.potential_energy()
    @property 
    def kinetic_energy(self):
        return self.gravity_code.kinetic_energy()
    def stop(self):
        if self.stellar_code != None:
            self.stellar_code.stop()
        self.gravity_code.stop()
        pass

    def read_perturber_list(self, filename):
        self.perturber_list = read_set_from_file(filename, close_file=True)
        self.perturber_list = list(self.perturber_list.iter_history())
        self.perturber_list_index = len(self.perturber_list)-1
        
        pstars = self.perturber_list[self.perturber_list_index-1].copy()
        dt_perturbation = pstars[pstars.name=="Sun"][0].age
        pstars = self.perturber_list[self.perturber_list_index].copy()
        star = pstars[pstars.name=="Sun"][0].copy()
        self.particles.add_particle(star)
        
        #print("N==", len(pstars), dt_perturbation.in_(units.Myr))
        if self.nperturbers<0:
            self.nperturbers = len(pstars)-1
        if self.nperturbers>len(pstars)-1:
            print("Too many nearest neighbors requested for cluster.")
            exit(-1)
            
    def add_perturbers(self):
        if self.perturber_list_index==0:
            print("No perturbers added because end time reached.")
            time_next = self.get_next_perturber_time()
            print(f"Integrate without perturbers: from time={self.model_time.value_in(units.Myr)} to {(time_next).value_in(units.Myr)}")
            return
        
        perturbers = self.perturber_list[self.perturber_list_index]
        star = perturbers[perturbers.name=="Sun"][0]
        #print("perturbed stars position=", star.x.in_(units.pc))
        #print("integrated position=", self.get_parent_star().x.in_(units.pc))
        perturbers -= star
        self.perturbers.add_particles(perturbers[:self.nperturbers])
        #print("Add perturbers:", len(self.perturbers))

        #print("before x=", self.perturbers.x.in_(units.pc))
        sun = self.particles[self.particles.name=="Sun"][0]
        self.particles.age = sun.age
        #self.particles.mass = sun.mass
        self.particles.position -= sun.position
        self.particles.velocity -= sun.velocity
        self.particles.position += star.position
        self.particles.velocity += star.velocity
        self.to_gravity.copy()

        self.gravity_code.particles.add_particles(self.perturbers)
        self.to_perturbers = self.gravity_code.particles.new_channel_to(self.perturbers)
        self.from_perturbers = self.perturbers.new_channel_to(self.gravity_code.particles)

        time_next = self.get_next_perturber_time()
        print(f"Integrate with perturbers: N={len(self.perturbers)} from time={self.model_time.value_in(units.Myr)} to {(time_next).value_in(units.Myr)}")

    def remove_perturbers(self):
        #print("After x=", self.perturbers.x.in_(units.pc))
        if len(self.perturbers)==0:
            print("No perturbers")
            return
        #print("Remove perturbers:", len(self.perturbers))
        self.gravity_code.particles.remove_particles(self.perturbers)
        self.perturbers.remove_particles(self.perturbers)
        #self.particles.remove_particles(self.perturbers)

    def remove_lost_planets(self):

        sun = self.particles[self.particles.name=="Sun"][0]
        panda = self.particles-sun
        asteroids = panda[panda.name=="asteroid"]
        planets = panda - asteroids
        get_orbital_elements_of_planetary_system(sun, planets)
        get_orbital_elements_of_planetary_system(sun, asteroids)

        lost_planets = panda[panda.eccentricity>1]
        la = Particles()
        if len(lost_planets)>0:
            p = lost_planets[lost_planets.type=="planet"]
            a = lost_planets[lost_planets.type=="asteroid"]
            print(f"At time={self.model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) unbound (planets, asteroids)")
            
            #print(lost_planets)
            lp = lost_planets.get_intersecting_subset_in(self.gravity_code.particles)
            la.add_particles(lp[lp.name=="asteroid"])
            self.gravity_code.particles.remove_particles(lp)
            self.particles.remove_particles(lost_planets)

        lost_planets = panda[np.isnan(panda.eccentricity)]
        if len(lost_planets)>0:
            p = lost_planets[lost_planets.type=="planet"]
            a = lost_planets[lost_planets.type=="asteroid"]
            print(f"At time={self.model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) nan (planets, asteroids)")
        
            #print(lost_planets)
            lp = lost_planets.get_intersecting_subset_in(self.gravity_code.particles)
            #lp = self.gravity_code.particles.get_intersecting_subset_in(lost_planets)
            #print(lp)
            la.add_particles(lp[lp.name=="asteroid"])
            self.gravity_code.particles.remove_particles(lp)
            self.particles.remove_particles(lost_planets)

        merged_planets = Particles()
        for pi in panda:
            #print("pi=", pi.eccentricity)
            if pi.semimajor_axis*(1-pi.eccentricity)<100|units.RSun:
                merged_planets.add_particle(pi)

        if len(merged_planets)>0:
            p = merged_planets[merged_planets.type=="planet"]
            a = merged_planets[merged_planets.type=="asteroid"]
            print(f"At time={self.model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) merged (planets, asteroids)")
        
            mp = merged_planets.get_intersecting_subset_in(self.gravity_code.particles)
            self.gravity_code.particles.remove_particles(mp)
            self.particles.remove_particles(merged_planets)

        print(f"Time={self.model_time.in_(units.Myr)}: a={planets.semimajor_axis.in_(units.au)}, e={planets.eccentricity}")


            
    def write_planetary_system(self):
        self.particles.age = self.model_time
        star = self.particles[self.particles.name=="Sun"][0]
        filename = "lp_planets_key_"+str(star.key)+".amuse"
        #planets = self.particles-star
        #print(f"Write file {filename}") 
        write_set_to_file(self.particles,
                          filename,
                          close_file=True,
                          append_to_file=True)
        print("File written:", filename, "N_planets=", len(self.particles))

            
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


def integrate_planetary_system(Nasteroids):

    galaxy_code = MilkyWay_galaxy()
    Rinit = 8.5|units.kpc
    
    star = Particle()
    star.type = "star"
    star.name = "Sun"
    star.mass = 1|units.MSun
    star.position = (0,0,0) | units.au
    star.velocity = (0,0,0) | units.kms
    
    star.x += Rinit
    star.vy = 1.0*galaxy_code.circular_velocity(Rinit)
    
    planetary_system = generate_planetary_system(star.as_set(),
                                                 Nasteroids)

    minor_bodies = planetary_system-star
    planets = minor_bodies[minor_bodies.type=="planet"]
    asteroids = minor_bodies-planets
    get_orbital_elements_of_planetary_system(star, planets)
    print(planets.semimajor_axis.in_(units.au))

    cluster_code = PlanetarySystemIntegrationWithPerturbers(planetary_system, Nnn)
    
    cluster_code.add_particles(planetary_system)
    
    gravity = bridge.Bridge(verbose=False)
    gravity.add_system(cluster_code, (galaxy_code,))
    gravity.timestep=0.01|units.Myr
    
    model_time = 0|units.Myr
    write_planetary_system(model_time, planetary_system)

    sys.stdout.flush()
    
    dt_diag = 1 |units.Myr
    t_diag = dt_diag
    t_end = 100|units.Myr
    dt = 0.01 | units.Myr
    while model_time<t_end:
        model_time += dt

        gravity.evolve_model(model_time)
        print(f"Time={model_time.in_(units.Myr)}: N= {len(gravity.particles)}")
        
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
    result.add_option("-N", dest="Nasteroids", type="int",default = 100,
                      help="number of asteroids [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    np.random.seed(31415)
    integrate_planetary_system(o.Nasteroids)
    
