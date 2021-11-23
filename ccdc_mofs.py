#!/usr/bin/env python
'''
https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/substructure_searching.html

Substructure searching with geometric measurements

MOF subset is here
/opt/ohpc/pub/ccdc/CSD_2021/csd/subsets/


'''
import time
import ccdc.search
import ccdc.io
from ccdc.io import CrystalReader, MoleculeWriter, CrystalWriter
from os.path import join
import os, sys
import csv
import numpy as np

# structure stuff
sys.path.insert(0, '/home/pboyd/modules/fa3ps')
from faps import Structure, Atom, Cell, minimum_image, min_distance

avo = 6.02214076e23 
cm3 = 1e24  # A^3 / cm^3

outdir = os.getcwd()
csd_subsetdir = join(os.environ['CSDHOME'], 'subsets')
datafile = 'csd_data.csv'
def adsorbaphore_search():
    searcher = ccdc.search.SubstructureSearch()
    benz_smarts = 'a1aaaaa1'
    arom1 = ccdc.search.SMARTSSubstructure(benz_smarts)
    arom2 = ccdc.search.SMARTSSubstructure(benz_smarts)
    
    sub1 = searcher.add_substructure(arom1)
    sub2 = searcher.add_substructure(arom2)
    
    
    searcher.add_plane('PLANE1', (sub1, 0), (sub1, 1), (sub1, 2),
                                 (sub1, 3), (sub1, 4), (sub1, 5))
    
    searcher.add_plane('PLANE2', (sub2, 0), (sub2, 1), (sub2, 2),
                                 (sub2, 3), (sub2, 4), (sub2, 5))
    
    searcher.add_centroid('CENT1', (sub1, 0), (sub1, 1), (sub1, 2),
                                   (sub1, 3), (sub1, 4), (sub1, 5))
    
    searcher.add_centroid('CENT2', (sub2, 0), (sub2, 1), (sub2, 2),
                                   (sub2, 3), (sub2, 4), (sub2, 5))
    
    # centroid between aromatic planes?
    searcher.add_centroid('CENT3', 'CENT1', 'CENT2')
    
    # Make sure the angle is between -5 and +5 deg?
    searcher.add_plane_angle_constraint('ANGLE', 'PLANE1', 'PLANE2', (-5, 5))
    searcher.add_distance_constraint('DIST', 'CENT1', 'CENT2', (6.8, 7.4), vdw_corrected=False, type='any')
    
    # Make sure the two aromatic planes are aligned
    searcher.add_vector('VEC1', 'CENT1', (sub1, 1))
    searcher.add_vector('VECN', 'CENT1', 'CENT2')
    searcher.add_vector_angle_constraint('ANGLE_N', 'VEC1', 'VECN', (81, 99)) 

    return searcher

searcher = adsorbaphore_search()
#searcher = benz_search()

# find hits
start_time = time.time()
# just get unique hits, compute number of hits per structure in the for loop below.

f = open(join(outdir, datafile), 'w')

cwriter = csv.writer(f)
cwriter.writerow(['CSD_NAME', 'UNIT_VOL_A^3', 'CRYSTAL_MOLAR_DENS_MMOL_CM^3', 'SUBSTRUCT_COUNT', 'SUBSTRUCT_DENS_MMOL_CM^3', 'PLANAR_ANGLE', 'PLANAR_ANGLE_STDEV', 'PLANAR_DIST', 'PLANAR_DIST_STDEV'])
# mine just the mofs
success_count, total_count = 0,0
searcher.settings.no_disorder = 'all'
searcher.settings.max_r_factor = 5.0
# mine just the mofs
mof_csd = join(csd_subsetdir, 'MOF_subset.gcd')
tot = len(ccdc.io.EntryReader(mof_csd))
for h in ccdc.io.EntryReader(mof_csd):
    # create a faps structure
    fstr = Structure(h.identifier)
    # need to ensure that h.molecule will give you the whole unit cell
    # atoms instead of just the asymmetric ones.
    for atom in h.molecule.atoms:
        if atom.coordinates is not None:
            coords = (atom.coordinates.x, atom.coordinates.y, atom.coordinates.z)
            fstr.atoms.append(Atom(idx=atom.index ,
                                   pos=coords, 
                                   at_type=atom.atomic_symbol, 
                                   mass=atom.atomic_weight,
                                   parent=fstr))
            assert (atom.index == len(fstr.atoms)-1)
        else:
            # just append the last known coordinates (prob a terrible idea :P)
            fstr.atoms.append(Atom(idx=atom.index,
                                   pos=coords, 
                                   at_type=atom.atomic_symbol, 
                                   mass=atom.atomic_weight,
                                   parent=fstr))

    cry = h.crystal
    pp = (cry.cell_lengths.a, cry.cell_lengths.b, cry.cell_lengths.c,
            cry.cell_angles.alpha, cry.cell_angles.beta, cry.cell_angles.gamma)
    fstr.cell.params = pp
    [at.get_fractional_coordinate() for at in fstr.atoms]
    local_hits = searcher.search(cry, max_hits_per_structure=10000)
    nhits = len(local_hits)
    #print('{0:s} hits: {1:d}'.format(h.identifier, nhits))
    if(nhits>0):
        rm = []

        # make sure nothing is in between the matched substructures
        # this does not work!
        # https://ccdc.cam.ac.uk/forum/csd_python_api/Crystallography/25cf2d8f-b11c-e711-84d4-005056975d8a#26cf2d8f-b11c-e711-84d4-005056975d8a
        # TODO(pboyd): make sure this works with a shifted benzene model (half in the unit cell + half out!)

        for i, hit in enumerate(local_hits):
            benz1, benz2 = hit.match_substructures()
            atom_list = hit.match_atoms()
            indices = [i.index for i in atom_list]
            # create first 
            c = [Atom(idx=None, pos=i.coordinates, at_type=i.atomic_symbol,mass=i.atomic_weight,parent=fstr) for i in benz1.atoms] 
            # create second
            c2 = [Atom(idx=None, pos=i.coordinates, at_type=i.atomic_symbol,mass=i.atomic_weight,parent=fstr) for i in benz2.atoms] 
            coords = ([i.coordinates[:] for i in benz1.atoms] + [i.coordinates[:] for i in benz2.atoms])
            # coords = [i.coordinates for i in atom_list] # atom list doesn't work! only in unit cell. must use substructures from match.

            # COM make sure periodic image shifted. for each aromatic ring
            centre = np.mean(coords, axis=0)

            # make sure middle of COMs is shifted.

            # measure minimum image distance of all non-adsorbaphore atoms
            # make sure not in van der waals radii.

            other_atoms = [fstr.atoms[i.index] for i in h.molecule.atoms if i.index not in indices]
            # i.fractional_coordinates, i.vdw_radius, i.coordinates
            # parent=fstr does this add the atom to fstr?
            comat = Atom(idx=None, pos=centre, at_type='X', mass=0.0, parent=fstr)
            # taking a faps Structure Atom, getting the vdw_radius from
            # hit.molecule.atom
            dists = [min_distance(comat, iat, cell=fstr.cell.cell)-
                    h.molecule.atoms[iat.idx].vdw_radius for iat in other_atoms]
            
            eval_ = [i > 0.0 for i in dists]

            # vector and find all atoms closest points to vector.
            # Benchmark this.

            # np.abs(np.cross(v, crystal_atoms_that_are_not_part_of_the_substructure - v[0]) / np.linalg.norm(v)) 
            # then subtract the vanderwaals radius of the crystal atoms.
            #sys.exit()
            if not np.all(eval_):
                nhits-=1
                rm.append(i)
        # delete matches that have stuff in between.
        rm.sort()
        rm.reverse()
        [local_hits.pop(i) for i in rm]
        # if there are any matches left, write to output.
    if(nhits>0):
        # spit out the matched atoms as a separate molecule file??
        # write_c2m_file breaks with error: AttributeError: 'NoneType' object has no attribute 'substructure_index'
        try:
            local_hits.write_c2m_file(join(outdir, '{0:s}_hits.c2m'.format(h.identifier)))
        except (AttributeError, NotImplementedError) as e:
            print('{0:s} yielded attribute error when writing c2m file.'.format(h.identifier))
        # count as a successful find
        success_count += 1
        a, d = zip(*[(s.constraints['ANGLE'], s.constraints['DIST']) for s in local_hits])
    
        AVG_ANGLE = np.absolute(np.array(a)).mean()
        STD_ANGLE = np.absolute(np.array(a)).std()
        AVG_DIST = np.array(d).mean()
        STD_DIST = np.array(d).std()
    
        vol = cm3 / cry.cell_volume * 1000 # how many unit cells fit into a cubic centimetre
        molar_dens = nhits * vol / avo # in mol/cm^3
        
        cwriter.writerow([h.identifier, cry.cell_volume, vol/avo, nhits, molar_dens, AVG_ANGLE, STD_ANGLE, AVG_DIST, STD_DIST])
        writer = ccdc.io.CrystalWriter(join(outdir, '{0:s}.cif'.format(h.identifier)), append=False)
        writer.write(h)
        writer.close()
    print('{0:s} hits: {1:d}'.format(h.identifier, nhits))
    total_count += 1
    print('{0:.2f}% complete ({1:d} of {2:d})'.format(total_count/tot*100., total_count, tot))
f.close()

end_time = time.time()
elapsed = end_time - start_time
f.close()
print ("Total found = {0:d}, Total MOFs = {1:d}, success rate = {2:.3f}".format(success_count, total_count, float(success_count)/float(total_count)))

