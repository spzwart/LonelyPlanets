from amuse.lab import *
from amuse.ext import solarsystem 

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

def make_new_planetary_system(Nasteroids):

    planetary_system = Particles()
    solar_system = new_solar_system()
    sun = planetary_system[planetary_system.name=="SUN"][0]
    sun.name = "Sun"
    solar_system.add_particles(sun)
    solar_system.position-=sun.position
    solar_system.velocity-=sun.velocity
    solar_system.remove_particle(sun)
    #planets = planetary_system[2:] #includes Jupiter and up
    planets = planetary_system[4:] # includes Earth and up
    planets = planets[:-1]
    planets.type = "planet"
    planetary_system.add_particles(planets)

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

    # rotate planetnary system
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

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Nasteroids", type="int",default = 100,
                      help="number of asteroids [%default]")
    result.add_option("-f", dest="outputfilename",
                      default = "initial_planetary_system.amuse",
                      help="output filename [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    np.random.seed(31415)
    system = make_new_planetary_system(o.Nasteroids)
    star = system[system.name=="Sun"]
    planets = system - star
    get_orbital_elements_of_planetary_system(star, planets)
    
    write_set_to_file(system, o.outfilename, close_file=True)    
