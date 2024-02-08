import numpy
import os.path
from amuse.lab import *

def new_rotation_matrix_from_euler_angles(phi, theta, chi):
    cosp=numpy.cos(phi)
    sinp=numpy.sin(phi)
    cost=numpy.cos(theta)
    sint=numpy.sin(theta)
    cosc=numpy.cos(chi)
    sinc=numpy.sin(chi)
    #see wikipedia: http://en.wikipedia.org/wiki/Rotation_matrix
    return numpy.array(
        [[cost*cosc, -cosp*sinc + sinp*sint*cosc, sinp*sinc + cosp*sint*cosc], 
         [cost*sinc, cosp*cosc + sinp*sint*sinc, -sinp*cosc + cosp*sint*sinc],
         [-sint,  sinp*cost,  cosp*cost]])



def rotate(position, velocity, phi, theta, psi): # theta and phi in radians
    Runit = position.unit
    Vunit = velocity.unit
    matrix = new_rotation_matrix_from_euler_angles(phi, theta, psi)
    return (numpy.dot(matrix, position.value_in(Runit)) | Runit,
           numpy.dot(matrix, velocity.value_in(Vunit)) | Vunit)

# select Eurler angles randomly. 

def random_Euler_angles():
    phi   = 2*numpy.pi*numpy.random.random()
    theta = numpy.arccos(1-2*numpy.random.random())
    chi   = 2*numpy.pi*numpy.random.random()
    return phi, theta, chi

def rotate_particle_set(bodies, phi, theta, chi):

    for bi in bodies:
        x_rot, v_rot = rotate(bi.position,
                              bi.velocity,
                              numpy.deg2rad(phi),
                              numpy.deg2rad(theta),
                              numpy.deg2rad(chi))
        bi.position = x_rot
        bi.velocity = v_rot
    return bodies

def test_rotate(phi, theta, chi):
    from amuse.ext.protodisk import ProtoPlanetaryDisk
    converter=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.au)
    disk = ProtoPlanetaryDisk(1000,
                              convert_nbody=converter,
                              Rmin=1,
                              Rmax=100,
                              q_out=1.0,
                              discfraction=0.01).result

    disk = rotate_particle_set(disk, phi, theta, chi)
    
    from matplotlib import pyplot

    pyplot.scatter(disk.x.value_in(units.au), disk.y.value_in(units.au))
    pyplot.xlabel("x-axis [au]")
    pyplot.ylabel("y-axis [au]")
    pyplot.show()
    pyplot.scatter(disk.x.value_in(units.au), disk.z.value_in(units.au))
    pyplot.xlabel("x-axis [au]")
    pyplot.ylabel("z-axis [au]")
    pyplot.show()

def rotate_bodies_isotropically(particles):
    phi, theta, chi = random_Euler_angles()
    phi = numpy.rad2deg(phi)
    theta = numpy.rad2deg(theta)
    chi = numpy.rad2deg(chi)
    disk = rotate_particle_set(particles, phi, theta, chi)
    return disk

def rotate_minor_bodies_around_star(star, bodies):
    disk = bodies[bodies.host==star.key]
    print(disk)
    disk.position -= star.position
    disk.velocity -= star.velocity
    disk = rotate_bodies_isotropically(disk)
    disk.position += star.position
    disk.velocity += star.velocity

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile", default = None,
                      help="output filename [%default]")
    result.add_option("-p",
                      dest="phi", type="float", default = None,
                      help="rotate under x-axis [%default]")
    result.add_option("-t",
                      dest="theta", type="float", default = -62.5,
                      help="rotate under y-axis [%default]")
    result.add_option("--type",
                      dest="rotate_type", default = "",
                      help="rotate particle type [%default]")
    result.add_option("-c",
                      dest="chi", type="float", default = 0.0,
                      help="rotate under z-axis [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    if o.filename==None:
        test_rotate(o.phi, o.theta, o.chi)
        exit

    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)
    if len(o.rotate_type)>0:
        stars = bodies[bodies.type==o.rotate_type]
        for star in stars:
            rotate_minor_bodies_around_star(star, bodies)
    else:
        com = bodies.center_of_mass()
        comv = bodies.center_of_mass_velocity()
        bodies.move_to_center()

        if o.phi==None:
            phi, theta, chi = random_Euler_angles()
            phi = numpy.rad2deg(phi)
            theta = numpy.rad2deg(theta)
            chi = numpy.rad2deg(chi)
            print("rotate randomly to (phi, theta, chi):", phi, theta, chi)
        else:
            phi = 0
            theta = 0
            chi = 0
        
        bodies = rotate_particle_set(bodies, phi, theta, chi)
        bodies.position += com
        bodies.velocity += comv

    if o.outfile == None:
        outfile = "particles_rotated.amuse"
    else:
        outfile = o.outfile
    index = 0
    time = 0|units.Myr
    write_set_to_file(bodies,
                      outfile, 'amuse',
                      timestamp=time,
                      overwrite_file=True,
                      version="2.0")
    

