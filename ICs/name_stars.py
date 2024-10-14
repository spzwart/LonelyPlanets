from amuse.lab import *

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("-F", 
                      dest="outfile",default = None,
                      help="output filename [%default]")
    result.add_option("--name", 
                      dest="name", 
                      default = "Sun",
                      help="disk mass [%default]")
    result.add_option("--nstars", 
                      dest="nstars", type="int",
                      default = 1,
                      help="number of stars with name [%default]")
    result.add_option("--seed", 
                      dest="seed", type="int",default = -1,
                      help="random number seed [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    if o.seed>0:
        numpy.random.seed(o.seed)
    else:
        print("random number seed from clock.")

    bodies = read_set_from_file(o.filename, 'hdf5', close_file=True)
    random_stars = bodies.random_sample(o.nstars)
    random_stars.name = o.name

    time = 0|units.Myr
    if o.outfile==None:
        filename = "added_names_to_subset_of_stars.amuse"
    else:
        filename = o.outfile
    write_set_to_file(bodies,
                      filename, 'amuse',
                      timestamp=time,
                      append_to_file=False,
                      version="2.0")
