"""
   Minimalistic routine for running a gravity code
"""
import os
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

#from matplotlib import pyplot as plt
import numpy as np
import math

class MilkyWay_galaxy(object):

    def __init__(self, Mb=1.40592e10| units.MSun,
                 Md=8.5608e10| units.MSun,
                 Mh=1.07068e11 | units.MSun):
        self.Mb = Mb
        self.Md = Md
        self.Mh = Mh

    def get_potential_at_point(self,eps,x,y,z):
        r = (x**2+y**2+z**2)**0.5
        R = (x**2+y**2)**0.5
        # bulge
        b1 = 0.3873 |units.kpc
        pot_bulge = -constants.G*self.Mb/(r**2+b1**2)**0.5 
        # disk
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        pot_disk = -constants.G*self.Md/(R**2+(a2+(z**2+b2**2)**0.5)**2)**0.5
        #halo
        a3 = 12.0 |units.kpc
        cut_off = 100 |units.kpc
        d1 =  r/a3
        c = 1+ (cut_off/a3)**1.02
        pot_halo = -constants.G*(self.Mh/a3)*d1**1.02/(1+ d1**1.02) \
                   - (constants.G*self.Mh/(1.02*a3))\
                     * (-1.02/c +numpy.log(c) + 1.02/(1+d1**1.02) \
                        - numpy.log(1.0+d1**1.02))
        return 2*(pot_bulge+pot_disk+ pot_halo) # multiply by 2 for
    						# a rigid potential

    def get_gravity_at_point(self, eps, x,y,z): 
        r = (x**2+y**2+z**2)**0.5
        R = (x**2+y**2)**0.5
        #bulge
        b1 = 0.3873 |units.kpc
        force_bulge = -constants.G*self.Mb/(r**2+b1**2)**1.5 
        #disk
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        d = a2+ (z**2+ b2**2)**0.5
        force_disk =-constants.G*self.Md/(R**2+ d**2 )**1.5
        #halo
        a3 = 12.0 |units.kpc
        d1 = r/a3
        force_halo = -constants.G*self.Mh*d1**0.02/(a3**2*(1+d1**1.02))
       
        ax = force_bulge*x + force_disk*x + force_halo*x/r
        ay = force_bulge*y + force_disk*y + force_halo*y/r
        az = force_bulge*z + force_disk*d*z/(z**2 + b2**2)**0.5 \
                           + force_halo*z/r 

        return ax,ay,az

    def circular_velocity(self,r):  
        z = 0 | units.kpc 
        b1 = 0.3873 |units.kpc
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        a3 = 12.0 |units.kpc

        rdphi_b = constants.G*self.Mb*r**2/(r**2+b1**2)**1.5
        rdphi_d = constants.G*self.Md*r**2/(r**2+(a2+(z**2+b2**2)**0.5)**2)**1.5
        rdphi_h = constants.G*self.Mh*(r/a3)**0.02*r/(a3**2*(1+(r/a3)**1.02))

        vel_circb = rdphi_b
        vel_circd = rdphi_d
        vel_circh = rdphi_h

        return (vel_circb + vel_circd + vel_circh)**0.5 

def resolve_supernova(supernova_detection, bodies, time):
    if supernova_detection.is_set():
        print("At time=", time.in_(units.Myr), \
              len(supernova_detection.particles(0)), 'supernova(e) detected')

        Nsn = 0
        for ci in range(len(supernova_detection.particles(0))):
            print(supernova_detection.particles(0))
            particles_in_supernova \
                = Particles(particles=supernova_detection.particles(0))
            natal_kick_x = particles_in_supernova.natal_kick_x
            natal_kick_y = particles_in_supernova.natal_kick_y
            natal_kick_z = particles_in_supernova.natal_kick_z
            print("Kick velocity:",
                  natal_kick_x.in_(units.kms),
                  natal_kick_y.in_(units.kms),
                  natal_kick_z.in_(units.kms))
                  
            particles_in_supernova \
                = particles_in_supernova.get_intersecting_subset_in(bodies)
            particles_in_supernova.vx += natal_kick_x
            particles_in_supernova.vy += natal_kick_y
            particles_in_supernova.vz += natal_kick_z
            Nsn += 1
        print('Resolved', Nsn, 'supernova(e)')
    
def ZAMS_radius(mass):
    log_mass = np.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun))**2
    alpha = 0.08353 + 0.0565*log_mass
    beta  = 0.01291 + 0.2226*log_mass
    gamma = 0.1151 + 0.06267*log_mass
    r_zams = pow(mass.value_in(units.MSun), 1.25) * (0.1148 + 0.8604*mass_sq) / (0.04651 + mass_sq)
    return r_zams | units.RSun

def merge_nn_with_sp(sp, nn, Nnn):
    perturbers = Particles()
    #perturbers.add_particles(nn)
    for spi, nni in zip(sp, nn):
        if nni not in perturbers:
            perturbers.add_particle(nni)
            if len(perturbers)>=Nnn:
                break
        if spi not in perturbers:
            perturbers.add_particle(spi)
            if len(perturbers)>=Nnn:
                break
    return perturbers

def find_perturbers(star, bodies, Nnn):
    nn = nearest_neightbors(bodies, star, Nnn)
    sp = strongest_perturber(bodies, star, Nnn)
    perturbers = merge_nn_with_sp(sp, nn, Nnn)
    return perturbers

def strongest_perturber(bodies, star, Nnn):
    r = ((bodies-star).position-star.position).lengths().value_in(units.pc)
    m = (bodies-star).mass.value_in(units.MSun)
    f = m/r**2
    f, p = (list(t) for t in zip(*sorted(zip(f, bodies-star))))
    #print(f, f[-Nnn:])
    b = Particles()
    for pi in p[-Nnn:]:
        b.add_particle(pi)
    return b

def nearest_neightbors(bodies, star, Nnn):
    
    r = ((bodies-star).position-star.position).lengths().value_in(units.pc)
    r, p = (list(t) for t in zip(*sorted(zip(r, bodies-star))))
    b = Particles()
    for pi in p[:Nnn]:
        b.add_particle(pi)
    #print("nn=", r, r[:Nnn], (b.position-star.position).lengths().in_(units.pc))
    return b

def get_filename(star):
    return "lps_key_"+str(star.key)+".amuse"

def write_perturbers_to_file(model_time, star, perturbers):
    star.age = model_time
    filename = get_filename(star)
    bodies = Particles()
    bodies.add_particle(star)
    bodies.add_particles(perturbers)
    write_set_to_file(bodies, filename,
                      close_file=True,
                      append_to_file=True)

def write_escaping_particles(escapers, model_time):
    escapers.age = model_time
    filename = "cl_escapers.amuse"
    write_set_to_file(escapers,
                      filename,
                      close_file=True,
                      append_to_file=True)

def find_escapers(bodies, model_time):
    rmax = 10|units.kpc
    escapers = bodies.select(lambda r: r.length()>rmax, ["position"])
    if (len(escapers))>0:
       write_escaping_particles(escapers, model_time)
    return escapers

def run_LonelyPlanets(bodies, 
                      time_end=10 | units.Myr,
                      dt=0.001|units.Myr,
                      Nnn=3):

    suns = bodies[bodies.name=="Sun"]
    print("N=", len(bodies), len(suns))

    galaxy_code = MilkyWay_galaxy()

    Rinit = 8500|units.pc
    bodies.x += Rinit
    bodies.vy += galaxy_code.circular_velocity(Rinit)
    bodies.vx += -10.1|units.kms
    bodies.vz += 7.5|units.kms

    stellar = SeBa()
    stellar.particles.add_particles(bodies)
    channel_from_se = stellar.particles.new_channel_to(bodies)
    channel_to_se = bodies.new_channel_to(stellar.particles)

    supernova_detection = stellar.stopping_conditions.supernova_detection
    supernova_detection.enable()
    
    stellar.evolve_model(0|units.Myr)
    channel_from_se.copy()
    index = 0
    fcluster = "cluster_i{0:06}.amuse".format(index)
    write_set_to_file(bodies, fcluster,
                      close_file=True,
                      append_to_file=True)
    
    converter=nbody_system.nbody_to_si(bodies.mass.sum(), time_end)
    
    cluster_code = ph4(converter, number_of_workers=12)
    #cluster_code = Petar(converter, number_of_workers=12)
    print(cluster_code.parameters)
    #cluster_code.parameters.epsilon_squared = (100|units.au)**2
    #cluster_code.parameters.timestep_parameter = 0.03
    cluster_code.particles.add_particles(bodies)
    energy_tot_init = cluster_code.kinetic_energy + cluster_code.potential_energy
    
    collision_detection = cluster_code.stopping_conditions.collision_detection
    collision_detection.enable()
    
    channel_from_gd = cluster_code.particles.new_channel_to(bodies)
    channel_to_gd = bodies.new_channel_to(cluster_code.particles)

    gravity = bridge.Bridge(verbose=False)
    gravity.add_system(cluster_code, (galaxy_code,))
    gravity.timestep=0.01|units.Myr

    model_time = 0 | units.Myr
    for star in suns:
        perturbers = find_perturbers(star, bodies, Nnn)
        write_perturbers_to_file(model_time, star, perturbers)

    t_diag = model_time
    dt_diag = 1|units.Myr
    Nsn = 0
    while gravity.model_time<time_end:

        stellar.evolve_model(model_time+dt)
        if supernova_detection.is_set():        
            resolve_supernova(supernova_detection, bodies, model_time)
            channel_from_se.copy_attributes(["mass"])
            channel_to_gd.copy_attributes(["mass", "vx", "vy", "vz"])
            Nsn += 1
            sn_file = "supernova_i{0:06}.amuse".format(Nsn)
            write_set_to_file(bodies, sn_file,
                              close_file=True,
                              append_to_file=True)
            
        else:
            channel_from_se.copy()
            channel_to_gd.copy_attributes(["mass"])
        model_time = stellar.model_time
            
        gravity.evolve_model(model_time)
        channel_from_gd.copy()

        escapers = find_escapers(bodies, model_time)
        if len(escapers)>0:
            print(f"At time={model_time.in_(units.Myr)}",
                  f" removed {len(escapers)} removed from Galaxy")
            bodies.remove_particles(escapers)
            bodies.synchronize_to(stellar.particles)
            bodies.synchronize_to(cluster_code.particles)
        
        for star in suns:
            perturbers = find_perturbers(star, bodies, Nnn)
            write_perturbers_to_file(model_time, star, perturbers)
        print(f"time={model_time.in_(units.Myr)}")

        if model_time>t_diag:
            t_diag += dt_diag
            index += 1
            fcluster = "cluster_i{0:06}.amuse".format(index)
            write_set_to_file(bodies, fcluster,
                              close_file=True,
                              append_to_file=True)

        if os.path.isfile("STOP"):
            print(f"Stop the simulation at time={model_time.in_(units.Myr)}")
            os.remove("STOP")
            break
            
    gravity.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Nnn", dest="Nnn", type="int",default = 4,
                      help="number of nearest neighbors [%default]")
    result.add_option("--NSuns", 
                      dest="NSuns", type="int",default = 1,
                      help="Number of Suns [%default]")
    result.add_option("--name", 
                      dest="name", default = "Sun",
                      help="identify Sun as planet-bearing star [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float",default = 1|units.Myr,
                      help="end time [%default]")
    result.add_option("--dt", unit=units.Myr,
                      dest="dt", type="float",default = 0.01|units.Myr,
                      help="time step [%default]")
    result.add_option("-f", 
                      dest="filename", default = "initial_cluster.amuse",
                      help="input filename [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    np.random.seed(31415)
    bodies = read_set_from_file(o.filename)
    if len(o.name)>0:
        suns = bodies[bodies.name==o.name]
    else:
        suns = bodies[bodies.mass>0.8|units.MSun]
        suns = suns[suns.mass<1.2|units.MSun]
        if len(suns)<o.NSuns:
            print("Too many suns requested in input snapshot.")
            exit(-1)
        if o.NSuns>0:
            suns = suns.random_sample(o.NSuns)

    print(np.mean(suns.mass.in_(units.MSun)).in_(units.MSun))
    suns.mass = 1|units.MSun
    suns.name = "Sun"
    print("N=", len(bodies))
    run_LonelyPlanets(bodies,
                      o.t_end,
                      o.dt,
                      Nnn=o.Nnn)

