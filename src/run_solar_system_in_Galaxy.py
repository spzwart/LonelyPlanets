"""
   Minimalistic routine for running a gravity code
"""
import sys
import time as wallclock
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

from run_cluster import MilkyWay_galaxy
from make_initial_cluster import ZAMS_radius
class ConstantParameters:
    __slots__ = ()
    MAXIMAL_TIMESTEP = 1 | units.Myr


class PlanetarySystemIntegrationWithPerturbers(object):
    def __init__(self,
                 maximal_timestep=1|units.Myr,
                 nperturbers = -1):
        self.model_time = 0 | units.yr
        self.restart_time = 0 | units.Myr
        self.gravity_code = None
        self.stellar_code = None
        self.collision_detection = None
        self.particles = Particles(0)
        self.nperturbers = nperturbers
        self.perturbers = Particles(0)
        self.perturber_list = None
        self.perturber_list_index = 0
        #self.collision_pairs = []

        self.key = None # host stellar key

        self.maximal_timestep = maximal_timestep
        self.minimal_number_of_planets = 1
        self.minimal_number_of_asteroids = 50


        self.constants = ConstantParameters()
        
        self.wct_stellar = 0 | units.s
        self.wct_gravity = 0 | units.s
        self.wct_initialization = 0 | units.s
        self.wct_total = 0 | units.s

        
        self.converter=nbody_system.nbody_to_si(1|units.MSun,
                                                1000|units.au)


    #@property
    #def current_time(self):
    #    return self.restart_time + self.model_time
    
    def print_wallclock_time(self):
        print(f"Timing: at t={self.model_time.in_(units.Myr)}:",
              f" initialization={self.wct_initialization.value_in(units.s)}", 
              f" stellar={self.wct_stellar.in_(units.s)}", 
              f" gravity = {self.wct_gravity.in_(units.s)}",
              f" Total computing = {(self.wct_stellar+self.wct_gravity).in_(units.s)}")
        
    def get_perturbed_particle(self):
        perturbed_particle = self.particles[self.particles.name=="Sun"]
        if len(perturbed_particle)==1:
            return perturbed_particle[0]
        else:
            print("No or too many perturbed particle(s)")
            print("Stop")
            exit(-1)

    def initialize_sun_at_galactic_position(self, R_init):
        
        galaxy_code = MilkyWay_galaxy()
        if self.key == None:
            star = Particle()
        else:
            star = Particle(key=int(self.key))
        star.type = "star"
        star.name = "Sun"
        star.mass = 1|units.MSun
        star.position = (0,0,0) | units.pc
        star.velocity = (0,0,0) | units.kms
        star.x += R_init
        vcirc = galaxy_code.circular_velocity(R_init)
        star.vy += vcirc
        return star
        
    def add_planetary_system(self,
                             planetary_system=None,
                             Nasteroids=0):
        if planetary_system==None:
            self.add_new_planetary_system(Nasteroids)
        else:
            self.add_existing_planetary_system(planetary_system)
            
    def add_existing_planetary_system(self, planetary_system):
        
        sun = planetary_system[planetary_system.name=="Sun"]
        star_from_perturber_list = self.particles[self.particles.name=="Sun"][0]

        # move planetary system to perturber-list star
        planetary_system.position -= sun.position
        planetary_system.velocity -= sun.velocity
        planetary_system.position -= star_from_perturber_list.position
        planetary_system.velocity -= star_from_perturber_list.velocity
        #And remove the perturber-list star
        self.particles.remove_particle(star_from_perturber_list)
        
        self.particles.add_particles(planetary_system)
        self.key = self.get_perturbed_particle().key

    def add_new_planetary_system(self, Nasteroids):

        parent_star = Particles()
        if len(self.particles)>0:
            parent_star.add_particle(self.get_perturbed_particle())
        else:
            Rinit = 8500 | units.pc
            sun = self.initialize_sun_at_galactic_position(Rinit)
            parent_star.add_particle(sun)
            self.particles.add_particles(parent_star)

        planetary_system = new_solar_system()
        sun = planetary_system[planetary_system.name=="SUN"][0]
        #parent_star.mass = sun.mass
        planetary_system.position-=sun.position
        planetary_system.velocity-=sun.velocity
        planetary_system.remove_particle(sun)
        #planets = planetary_system[2:] #includes Jupiter and up
        planets = planetary_system[4:] # includes Earth and up
        planets = planets[:-1]
        planets.type = "planet"

        if self.minimal_number_of_planets<=len(planets):
            self.minimal_number_of_planets = len(planets)-self.minimal_number_of_planets
        converter=nbody_system.nbody_to_si(parent_star.mass.sum(), 1|units.au)
        asteroids = ProtoPlanetaryDisk(Nasteroids,
                                       densitypower=1.5, Rmin=1, Rmax=100,
                                       q_out=1, discfraction=0.01,
                                       convert_nbody=converter).result
        asteroids.mass = 0 | units.MEarth
        asteroids.name = "asteroid"
        asteroids.type = "asteroid"

        # rotate planetary system
        from rotate import rotate_particle_set
        phi = 0 #| units.deg rotate over x
        theta = -60.2 #| units.deg  rotate over y
        chi = 0 #| units.deg rotate around z
        rotate_particle_set(planets, phi, theta, chi)
        rotate_particle_set(asteroids, phi, theta, chi)
        print("rotate planetary system along the y-axis")
        
        # translate planetary system
        planets.position += parent_star.position
        planets.velocity += parent_star.velocity
        asteroids.position += parent_star.position
        asteroids.velocity += parent_star.velocity
    
        self.particles.add_particles(planets)
        self.particles.add_particles(asteroids)

        self.key = self.get_perturbed_particle().key

    def determine_orbital_parameters(self):
        #parent_star = self.particles[self.particles.type=="star"]
        parent_star = self.particles[self.particles.name=="Sun"]
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
    def start_gravity_code(self, integrator=0):
        self.gravity_code = Huayno(convert_nbody=self.converter,
                                   redirection="none")
        self.gravity_code.parameters.inttype_parameter = integrator

        sun = self.particles[self.particles.name=="Sun"]
        #print("Add planetary (start_gravity_code) system:", sun)

        self.gravity_code.particles.add_particles(self.particles)
        #self.gravity_code.parameters.timestep_parameter = 0.015
        print(self.gravity_code.parameters)

        self.from_gravity = self.gravity_code.particles.new_channel_to(self.particles)
        self.to_gravity = self.particles.new_channel_to(self.gravity_code.particles)

        self.collision_detection = self.gravity_code.stopping_conditions.collision_detection
        self.collision_detection.enable()

    def get_last_perturber_time(self):
        if self.perturber_list==None:
            return 0 | units.Myr
        else:
            return self.perturber_list[0][0].age
    def get_next_perturber_time(self):
        print("Perturber list index:", self.perturber_list_index)
        if self.perturber_list_index>=1:
            t_next_perturber = self.perturber_list[self.perturber_list_index-1][0].age
        else:
            t_next_perturber = self.model_time + self.maximal_timestep
        print("Next perturber time=", t_next_perturber.in_(units.Myr))
        return t_next_perturber
    
    def get_current_perturber_time(self):
        return self.perturber_list[self.perturber_list_index][0].age

    def evolve_model(self, model_time):
        print("Evolve to", model_time.in_(units.Myr))

        if self.perturber_list != None:
            self.model_time = self.perturber_list[self.perturber_list_index][0].age
        print("Current model time=", self.model_time.in_(units.Myr))
        print("Evolve to time=", model_time.in_(units.Myr))
        while self.model_time<model_time:

            #print("time=", self.model_time.in_(units.Myr))

            if self.model_time<self.get_last_perturber_time():
                self.add_perturbers()

            self.print_status()

            self.evolve_system_for_one_step()

            self.remove_perturbers()
            
            self.remove_lost_planets()

    def evolve_system_for_one_step(self):
        print("Evolve from", self.model_time.in_(units.Myr),
              "to ", self.get_next_perturber_time().in_(units.Myr),
              "time index=", self.perturber_list_index)
        sun = self.particles[self.particles.name=="Sun"][0]
        print("pos=", self.model_time.in_(units.Myr), sun.position.in_(units.pc))
        if self.stellar_code != None:
            print("Evolve stars.")
            wct = wallclock.time() 
            self.stellar_code.evolve_model(self.model_time)
            self.wct_stellar += (wallclock.time()-wct) | units.s
            self.from_stellar.copy()

        # To test the collosion handling
        #blow up stellar radius to test collisions
        #sun = self.particles[self.particles.name=="Sun"]
        #sun.radius = 4|units.au
        #self.particles.radius = 1|units.au
        
        self.to_gravity.copy()
        #print("R=", self.gravity_code.particles.radius.in_(units.au))
        #Subtract the restart time, because gravity cannot start
        #at any random moment in time
        time_next = self.get_next_perturber_time()-self.restart_time
        while self.gravity_code.model_time<time_next-(1|units.yr):
            #print("time=", self.gravity_code.model_time.in_(units.Myr)-time_next.in_(units.Myr), self.gravity_code.model_time.in_(units.Myr)-self.model_time.in_(units.Myr))
            wct = wallclock.time()
            self.gravity_code.evolve_model(time_next)
            self.wct_gravity += (wallclock.time()-wct) | units.s
            self.flag_colliding_stars()
            self.from_gravity.copy()
        
        if len(self.perturbers)>0:
            self.to_perturbers.copy()
        if self.perturber_list_index==0:
            print("End of perturber list reached at time=",
                  self.model_time.in_(units.Myr))
        else:
            self.perturber_list_index -= 1

        self.model_time = time_next + self.restart_time
        #self.get_current_perturber_time()
        #print("Evolution done, current index=", self.perturber_list_index)
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
    def terminating_criterium_reached(self):
        planets = self.particles[self.particles.type=="planet"]
        asteroids = self.particles[self.particles.type=="asteroid"]
        stop_run = False
        if len(planets)<=self.minimal_number_of_planets:
            stop_run = True
        if len(asteroids)<=self.minimal_number_of_asteroids:
            stop_run = True
        return stop_run
    
    def stop(self):
        if self.stellar_code != None:
            self.stellar_code.stop()
        self.gravity_code.stop()
        pass

    def set_perturber_list_index(self):
        self.perturber_list_index = len(self.perturber_list)-1
        print("Perturber list time:", self.model_time.in_(units.Myr))
        print("Perturber list index:", self.perturber_list_index)
        if self.model_time>0|units.Myr:
            for pli in range(self.perturber_list_index):
                #print("Perturber list index at age=",
                #      self.perturber_list[pli][0].age.in_(units.Myr),
                #      self.perturber_list[pli][0].position.in_(units.pc))
                if self.perturber_list[pli][0].age<=self.model_time:
                    self.perturber_list_index = pli
                    break
            print("Start calculation from perturber list time=", self.perturber_list[pli][0].age.in_(units.Myr), "index=", self.perturber_list_index)
    
    def read_perturber_list(self, filename):
        if self.nperturbers==0:
            key = filename.split("_")[-1].split(".amuse")[0]
            self.set_stellar_identity(key)
            #if self.nperturbers<0 or self.nperturbers>=1:
        else:
            self.perturber_list = read_set_from_file(filename,
                                                     close_file=True)
            self.perturber_list = list(self.perturber_list.iter_history())

            self.set_perturber_list_index()
            
            pstars = self.perturber_list[self.perturber_list_index-1].copy()
            dt_perturbation = pstars[pstars.name=="Sun"][0].age
            pstars = self.perturber_list[self.perturber_list_index].copy()
            if self.nperturbers<0 or self.nperturbers>len(pstars):
                self.nperturbers = len(pstars)-1
            print(f"Number of perturbers adopted: N={self.nperturbers}")

            star = pstars[pstars.name=="Sun"][0].copy()
            self.set_stellar_identity(star.key)
            print("Star from perturber list: pos=", star.position.in_(units.pc))

            ## Add star from perturber list in particle set
            ## We should keep its position and velocity, 
            ## but replace it with the planetary system's host
            ## in add_existing_planetary_system
            self.particles.add_particle(star)
            
            print("Current particles:", self.particles)
            star_tmp = self.particles[self.particles.name=="Sun"][0]
            print("Star from planetary system: pos=", star_tmp.position.in_(units.pc))
        
            #print("N==", len(pstars), dt_perturbation.in_(units.Myr))
            if self.nperturbers<0:
                self.nperturbers = len(pstars)-1
            if self.nperturbers>len(pstars)-1:
                print("Too many nearest neighbors requested for cluster.")
                exit(-1)
            
    def add_perturbers(self):
        self.add_perturbers_self_centered()
        #self.add_perturbers_sun_centered()
        
    def add_perturbers_sun_centered(self):
        if self.perturber_list_index==0:
            print("No perturbers added because end time reached.")
            time_next = self.get_next_perturber_time()
            print(f"Integrate without perturbers: from time={self.model_time.value_in(units.Myr)}.")
            return
        
        perturbers = self.perturber_list[self.perturber_list_index]
        star = perturbers[perturbers.name=="Sun"][0]
        #print("perturbed stars position=", star.x.in_(units.pc))
        #print("integrated position=", self.get_parent_star().x.in_(units.pc))
        perturbers -= star
        perturbers.position -= star.position
        perturbers.velocity -= star.velocity

        sun = self.particles[self.particles.name=="Sun"][0]
        perturbers.position += sun.position
        perturbers.velocity += sun.velocity
        
        dmin = (perturbers.position-sun.position).lengths().min()
        r = (perturbers.position-sun.position).lengths()
        m = perturbers.mass
        fmax = (m/r**2).max()
        print("minium distance=", dmin.in_(units.pc), "fmax=", fmax)

        dlim = 1|units.pc
        flim = ((1|units.MSun)/dlim**2).max()
        print("limits in distance:", dmin/dlim, "and force:", fmax/flim)
        
        if dmin/dlim<1.0 or fmax/flim>10:
            self.perturbers.add_particles(perturbers[:self.nperturbers])
            print("Add perturbers:", len(self.perturbers))
            print("self.perturbers.position.in_(units.pc)")

            #print("before x=", self.perturbers.x.in_(units.pc))
            self.particles.age = sun.age
            #self.particles.mass = sun.mass
            #self.particles.position -= sun.position
            #self.particles.velocity -= sun.velocity
            #self.particles.position += star.position
            #self.particles.velocity += star.velocity
            self.to_gravity.copy()

            self.gravity_code.particles.add_particles(self.perturbers)
            self.to_perturbers = self.gravity_code.particles.new_channel_to(self.perturbers)
            self.from_perturbers = self.perturbers.new_channel_to(self.gravity_code.particles)

            time_next = self.get_next_perturber_time()
            print(f"Integrate with perturbers: N={len(self.perturbers)} from time={self.model_time.value_in(units.Myr)} to {(time_next).value_in(units.Myr)}")

    def add_perturbers_self_centered(self):
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
        sun = self.particles[self.particles.name=="Sun"][0]
        if self.perturber_list_index>=1:
            #print("before x=", self.perturbers.x.in_(units.pc))
            #self.particles.mass = sun.mass
            self.particles.position -= sun.position
            self.particles.velocity -= sun.velocity
            self.particles.position += star.position
            self.particles.velocity += star.velocity
            self.to_gravity.copy()
            
        self.particles.age = sun.age

        dmin = (perturbers.position-star.position).lengths().min()
        r = (perturbers.position-star.position).lengths()
        m = perturbers.mass
        fmax = (m/r**2).max()
        print("minium distance=", dmin.in_(units.pc),
              "fmax=", fmax.in_(units.kg/units.km**2),
              "mpert=", m.in_(units.MSun))

        dlim = 3.0|units.pc # rRHill of the Sun in Galactic potential
        flim = (0.1|units.MSun)/dlim**2
        #print("Limits in distance:", dmin/dlim, "and force:", fmax/flim)
        if dmin<dlim: 
            print(f"NN in perturbing distance ({dmin/dlim}), ", end='')
            if fmax>flim:            
                print(f"and strongly perturbing ({fmax/flim})")
            else:
                print(".")
        else:
            print(f"NN not in perturbing distance ({dmin/dlim}), ", end='')
            if fmax>flim:            
                print(f"but strongly perturbing ({fmax/flim})")
            else:
                print(".")
        
        if dmin<dlim or fmax>flim:
            self.perturbers.add_particles(perturbers[:self.nperturbers])

            self.gravity_code.particles.add_particles(self.perturbers)
            self.to_perturbers = self.gravity_code.particles.new_channel_to(self.perturbers)
            self.from_perturbers = self.perturbers.new_channel_to(self.gravity_code.particles)

        time_next = self.get_next_perturber_time()
        print(f"Integrate with perturbers: N={len(self.perturbers)} from time={self.model_time.value_in(units.Myr)} to {(time_next).value_in(units.Myr)}: dmin={dmin.in_(units.pc)}, fmax={fmax.in_(units.kg/units.km**2)}")



            
    def remove_perturbers(self):
        #print("After x=", self.perturbers.x.in_(units.pc))
        if len(self.perturbers)==0:
            print("No perturbers")
            return
        #print("Remove perturbers:", len(self.perturbers))
        self.gravity_code.particles.remove_particles(self.perturbers)
        self.perturbers.remove_particles(self.perturbers)
        #self.particles.remove_particles(self.perturbers)
        #self.particles.move_to_center()

        if self.perturber_list_index>=1:
            #set back on the Sun
            sun = self.particles[self.particles.name=="Sun"][0]
            self.particles.position -= sun.position
            self.particles.velocity -= sun.velocity
            #set back to the new position of the Sun
            new_perturbers = self.perturber_list[self.perturber_list_index-1]
            star = new_perturbers[new_perturbers.name=="Sun"][0]
            self.particles.position += star.position
            self.particles.velocity += star.velocity
        else:
            print("Simulation has run out of perturbers.")
            
        self.to_gravity.copy()        

    def remove_lost_planets(self):

        escapers = Particles()
        sun = self.particles[self.particles.name=="Sun"][0]
        panda = self.particles-sun
        asteroids = panda[panda.type=="asteroid"]
        planets = panda - asteroids
        get_orbital_elements_of_planetary_system(sun, planets)
        get_orbital_elements_of_planetary_system(sun, asteroids)

        lost_planets = panda[panda.eccentricity>1]
        la = Particles()
        if len(lost_planets)>0:
            lost_planets.name = "escaper"
            escapers.add_particles(lost_planets)
            p = lost_planets[lost_planets.type=="planet"]
            a = lost_planets[lost_planets.type=="asteroid"]
            print(f"At time={self.model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) unbound (planets, asteroids)")
            
            #print(lost_planets)
            #lp = lost_planets.get_intersecting_subset_in(self.gravity_code.particles)
            #la.add_particles(lp[lp.type=="asteroid"])
            self.gravity_code.particles.remove_particles(lost_planets)
            self.particles.remove_particles(lost_planets)

        lost_planets = panda[np.isnan(panda.eccentricity)]
        if len(lost_planets)>0:
            lost_planets.name = "unknown"
            escapers.add_particles(lost_planets)
            p = lost_planets[lost_planets.type=="planet"]
            a = lost_planets[lost_planets.type=="asteroid"]
            print(f"At time={self.model_time.in_(units.Myr)}: remove ({len(p)}, {len(a)}) nan (planets, asteroids)")
        
            #lp = lost_planets.get_intersecting_subset_in(self.gravity_code.particles)
            #la.add_particles(lp[lp.type=="asteroid"])
            self.gravity_code.particles.remove_particles(lost_planets)
            self.particles.remove_particles(lost_planets)


        # small pericenter
        pericentra = panda.semimajor_axis*(1-panda.eccentricity)
        mask = pericentra<10|units.RSun
        lost_particles = panda[mask]
        
        # near the Sun
        d = (sun.position-lost_particles.position).lengths()
        mask = d<0.1|units.au
        colliding_particles = lost_particles[mask]
        #print("Particles also close to the Sun: N=", len(lost_particles))
        print("Particles colliding with the Sun: N=",
              len(colliding_particles))

        # near pericenter
        if len(colliding_particles)>0:
            colliding_particles.name = "collision"
            escapers.add_particles(colliding_particles)
            p = colliding_particles[colliding_particles.type=="planet"]
            a = colliding_particles[colliding_particles.type=="asteroid"]
            print(f"At time={self.model_time.in_(units.Myr)}:",
                  f"remove ({len(p)}, {len(a)}) colliding (planets, asteroids)")
            lp = colliding_particles.get_intersecting_subset_in(self.gravity_code.particles)
            self.gravity_code.particles.remove_particles(lp)
            self.particles.remove_particles(colliding_particles)

        if len(escapers)>0:
            print(f"Number of particles lost N=: {len(escapers)}")
            self.write_escaping_particles(escapers)

        print(f"Time={self.model_time.in_(units.Myr)} Planet orbits: a={planets.semimajor_axis.in_(units.au)}, e={planets.eccentricity}")

        #print(f"Time={self.model_time.in_(units.Myr)} Asteroid orbits: a={np.sort(asteroids.semimajor_axis.value_in(units.au))}, e={np.sort(asteroids.eccentricity)}")

    def set_stellar_identity(self, key=None):
        self.key = key
        
    def get_stellar_identity(self):
        key = self.key
        if key==None:
            try:
                star = self.get_perturbed_particle()
                key = star.key
            except:
                print("No key found in perturbed particles.")
                key = "NoKey"
        return key
        
    def write_escaping_particles(self, escapers):
        escapers.age = self.model_time
        key = self.get_stellar_identity()
        filename = "lp_escapers_key_"+str(key)+".amuse"
        write_set_to_file(escapers,
                          filename,
                          close_file=True,
                          append_to_file=True)
        
    def write_planetary_system(self):
        self.particles.age = self.model_time
        key = self.get_stellar_identity()
        filename = "lp_planets_key_"+str(key)+".amuse"
        #planets = self.particles-star
        #print(f"Write file {filename}") 
        write_set_to_file(self.particles,
                          filename,
                          close_file=True,
                          append_to_file=True)
        print("File written:", filename, "N_planets=", len(self.particles))

    def merge_two_particles(self, particles_in_encounter):
        mass = particles_in_encounter.total_mass()
        com_pos = particles_in_encounter.center_of_mass()
        com_vel = particles_in_encounter.center_of_mass_velocity()
        if particles_in_encounter[0].mass>particles_in_encounter[1].mass:
            target = particles_in_encounter[0]
            projectile = particles_in_encounter[1]
        else:
            target = particles_in_encounter[1]
            projectile = particles_in_encounter[0]
        print(f"Merge {target.name} with {projectile.name} with masses={particles_in_encounter.mass.in_(units.MJupiter)}")
        
        target.mass = mass
        target.position = com_pos
        target.velocity = com_vel

        if "star" in target.type:
            target.radius = ZAMS_radius(target.mass)
            print(f"Collision with {target.name}.")
        elif "planet" in target.type:
            rho = (1|units.MJupiter)/(1|units.RJupiter)**3 
            target.radius = (target.mass/rho)**(1./3.)
            print("Collision with {target.name}.")
        else:
            print("Collision between two asteroids.")
        #self.particles.add_particles(new_particle)
        #self.particles.remove_particles(projectile.as_set())
        return projectile.as_set()

    def flag_colliding_stars(self):
        if self.collision_detection.is_set():
            E_coll = self.gravity_code.kinetic_energy + \
                self.gravity_code.potential_energy
            print("At time=", self.gravity_code.model_time.in_(units.Myr),
                  "number of encounters=",
                  len(self.collision_detection.particles(0)))
            Nenc = 0
            for ci in range(len(self.collision_detection.particles(0))):
                print("Collision between: ",
                      self.collision_detection.particles(0)[ci].key,
                      "and",
                      self.collision_detection.particles(1)[ci].key)
                particles_in_encounter \
                    = Particles(particles=[self.collision_detection.particles(0)[ci],
                                           self.collision_detection.particles(1)[ci]])
                particles_in_encounter \
                    = particles_in_encounter.get_intersecting_subset_in(self.particles)

                #collision_pairs.add_particles(particles_in_encounter)

                projectile = self.merge_two_particles(particles_in_encounter)
                # you cannot use synchronize_to here, because
                # the local particle set does not contain the perturbers
                # They would be removed from the gravity code
                #self.particles.synchronize_to(self.gravity_code.particles)
                self.particles.remove_particles(projectile)
                self.gravity_code.particles.remove_particles(projectile)
                self.to_gravity.copy()

                Nenc += 1
                print("Resolve encounter Number:", Nenc)
            dE_coll = E_coll - (self.gravity_code.kinetic_energy + self.gravity_code.potential_energy)
            print("dE_coll =", dE_coll, "N_enc=", Nenc)
        sys.stdout.flush()
        
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
    planets.inclination = inc
    planets = planets[planets.eccentricity>=1]
    planets.eccentricity = 10
    planets.semimajor_axis = 1.e+10 | units.au
    planets = planets[np.isnan(planets.eccentricity)]
    planets.eccentricity = 10
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

    wct_initialization = wallclock.time()
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

    cluster_code.wct_initialization = (wallclock.time()-wct_initialization) | units.s

    model_time = 0|units.Myr
    write_planetary_system(model_time, planetary_system)

    sys.stdout.flush()
    
    dt_diag = 0.1 |units.Myr
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
    cluster_code.print_wallclock_time()
    
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

    converter=nbody_system.nbody_to_si(parent_star.mass.sum(), 1|units.au)
    asteroids = ProtoPlanetaryDisk(Nasteroids,
                                   densitypower=1.5, Rmin=1, Rmax=10000,
                                   q_out=1, discfraction=0.01,
                                   convert_nbody=converter).result
    asteroids.mass = 0 | units.MEarth
    asteroids.name = "asteroid"
    asteroids.type = "asteroid"

    # rotate planetary system
    from rotate import rotate_particle_set
    phi = 90 | units.deg
    #theta = 60.2 | units.deg
    theta = 0 | units.deg
    chi = 0 | units.deg
    rotate_particle_set(planets, phi, theta, chi)
    rotate_particle_set(asteroids, phi, theta, chi)

    #translate planetary system
    planets.position += parent_star.position
    planets.velocity += parent_star.velocity
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
    
