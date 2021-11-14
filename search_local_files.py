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

avo = 6.02214076e23 
cm3 = 1e24  # A^3 / cm^3

outdir = os.getcwd()
struct_dir = sys.argv[1] 
csd_subsetdir = join(os.environ['CSDHOME'], 'subsets')
datafile = 'substructure_data.csv'

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

# make sure the angle is between -5 and +5 deg?
searcher.add_plane_angle_constraint('ANGLE', 'PLANE1', 'PLANE2', (0, 10))
searcher.add_distance_constraint('DIST', 'CENT1', 'CENT2', (6.5, 7.2), vdw_corrected=False, type='any')

# Make sure the two aromatic planes are aligned
searcher.add_vector('VEC1', 'CENT1', (sub1, 0))
searcher.add_vector('VECN', 'CENT1', 'CENT2')
searcher.add_vector_angle_constraint('ANGLE_N', 'VEC1', 'VECN', (80, 100)) 

# find hits
start_time = time.time()
# just get unique hits, compute number of hits per structure in the for loop below.

f = open(join(outdir, datafile), 'w')

cwriter = csv.writer(f)
cwriter.writerow(['CSD_NAME', 'UNIT_VOL_A^3', 'CRYSTAL_MOLAR_DENS_MMOL_CM^3', 'SUBSTRUCT_COUNT', 'SUBSTRUCT_DENS_MMOL_CM^3', 
                  'PLANAR_ANGLE', 'PLANAR_ANGLE_STDEV', 'PLANAR_DIST', 'PLANAR_DIST_STDEV'])
#searcher.settings.no_disorder = 'all'
#searcher.settings.max_r_factor = 5.0
# mine just the mofs
success_count, total_count = 0,0
files = [i for i in os.listdir(struct_dir) if i.endswith('.cif')]#[:30]

reader = CrystalReader([join(struct_dir, j) for j in files])

for h in reader:
    h.assign_bonds()
    local_hits = searcher.search(h, max_hits_per_structure=10000)
    nhits = len(local_hits)
    #print('{0:s} hits: {1:d}'.format(h.identifier, nhits))
    if(nhits>0):
        rm = []

        # make sure nothing is in between the matched substructures
        for i, hit in enumerate(local_hits):
            benz1, benz2 = hit.match_substructures()
            atom_list = hit.match_atoms()
            # evaluate if an atoms from the first molecule are in line of sight with the second
            eval_ = [a.is_in_line_of_sight(b) for a in atom_list[:6] for b in atom_list[6:]]
            if not np.all(eval_):
                nhits-=1
                rm.append(i)
        # delete matches that have stuff in between.
        rm.sort()
        rm.reverse()
        [local_hits.pop(i) for i in rm]
        # if there are any matches left, write to output.
        if(len(local_hits)>0):
            # spit out the matched atoms as a separate molecule file??
            # outfile = join(outdir, '{0:s}_hits.mol2'.format(h.identifier))
            # molwriter = ccdc.io.MoleculeWriter(outfile)
            # mols = local_hits.superimpose()
            # [molwriter.write(mm) for mm in mols]
            # molwriter.close()
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

            #writer = ccdc.io.CrystalWriter(join(outdir, '{0:s}.cif'.format(h.identifier)), append=False)
            #writer.write(h.crystal)
            vol = cm3 / h.cell_volume * 1000 # how many unit cells fit into a cubic centimetre
            molar_dens = nhits * vol / avo # in mol/cm^3
            
            cwriter.writerow([h.identifier, h.cell_volume, vol/avo, nhits, molar_dens, AVG_ANGLE, STD_ANGLE, AVG_DIST, STD_DIST])
    total_count += 1

end_time = time.time()
elapsed = end_time - start_time
f.close()
print ("Total found = {0:d}, Total MOFs = {1:d}, success rate = {2:.3f}".format(success_count, total_count, float(success_count)/float(total_count)))

