# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 10:30:35 2019

@author: 1203087
"""

import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.lib.distances

import numpy as np
#from itertools import izip

try:
    import matplotlib

    matplotlib.use('agg')  # no interactive plotting, only save figures
    import pylab

    have_matplotlib = True
except ImportError:
    have_matplotlib = False

#MDAnalysis.start_logging()
file='indentation/indentation/dump/aver50.15.data'
simulation=MDAnalysis.Universe(file, atom_style='id type x y z')
all_atoms=simulation.atoms
natom=all_atoms.n_atoms
#print(natom)

posits=all_atoms.positions
#print(posits.ndim, posits.shape, posits.dtype)
#print(posits[0,0],posits[0,1],posits[0,2])
#print(posits[250046,0],posits[250046,1],posits[250046,2])


# set up rdf
dmin, dmax = 0.0, 16.0
nbins = 100
rdf, edges = np.histogram([0], bins=nbins, range=(dmin, dmax))
rdf *= 0
rdf = rdf.astype(np.float64)  # avoid possible problems with '/' later on

for i in range(1,natom):    #loop from 0 ~ n-1
    selector= "bynum " + str(i)
    atom_i=all_atoms.select_atoms(selector)
    #print( "center atom is ", atom_i )
    selector= "bynum " + str(i+1) + ":" + str(natom)
    rest_atoms=all_atoms.select_atoms(selector)
    #print( "rest atom are ", rest_atoms )
    dist = np.zeros( (1, rest_atoms.n_atoms) , dtype=np.float64)
    MDAnalysis.lib.distances.distance_array( atom_i.positions, rest_atoms.positions, result=dist, backend='serial')
    new_rdf, edges = np.histogram(dist, bins=nbins, range=(dmin, dmax))
    rdf += new_rdf

boxvolume = 1
numframes = 1
boxvolume /= numframes  # average volume

# Normalize RDF
radii = 0.5 * (edges[1:] + edges[:-1])
vol = (4. / 3.) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))

# normalization to the average density n/boxvolume in the simulation
density = natom / boxvolume
# This is inaccurate when solutes take up substantial amount
# of space. In this case you might want to use
## import MDAnalysis.units
## density = MDAnalysis.units.convert(1.0, 'water', 'Angstrom^{-3}')
norm = density * (natom - 1) / 2 * numframes
rdf /= norm * vol

outfile = './rdf.dat'
with open(outfile, 'w') as output:
    for radius, gofr in zip(radii, rdf):
        output.write("{radius:8.3f} \t {gofr:8.3f}\n".format(**vars()))
#print( "g(r) data written to {outfile!r}" % (**vars()) )

if have_matplotlib:
    matplotlib.rc('font', size=14)
    matplotlib.rc('figure', figsize=(5, 4))
    pylab.clf()
    pylab.plot(radii, rdf, linewidth=3)
    pylab.xlabel(r"distance $r$ in $\AA$")
    pylab.ylabel(r"radial distribution function $g(r)$")
    pylab.savefig("./rdf.png")    