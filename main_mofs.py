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
from os.path import join
import os
import csv
import numpy as np

avo = 6.02214076e23 
cm3 = 1e24  # A^3 / cm^3

outdir = os.getcwd()
csd_subsetdir = join(os.environ['CSDHOME'], 'subsets')
datafile = 'csd_data.csv'

searcher = ccdc.search.SubstructureSearch()
benz_smarts = 'c1ccccc1'
arom1 = ccdc.search.SMARTSSubstructure(benz_smarts)
arom2 = ccdc.search.SMARTSSubstructure(benz_smarts)

sub1 = searcher.add_substructure(arom1)
sub2 = searcher.add_substructure(arom2)

searcher.add_plane('PLANE1', (sub1, 0), (sub1, 1), (sub1, 2))
searcher.add_plane('PLANE2', (sub2, 0), (sub2, 1), (sub2, 2))

searcher.add_centroid('CENT1', (sub1, 0), (sub1, 1), (sub1, 2),
                               (sub1, 3), (sub1, 4), (sub1, 5))

searcher.add_centroid('CENT2', (sub2, 0), (sub2, 1), (sub2, 2),
                               (sub2, 3), (sub2, 4), (sub2, 5))

# centroid between aromatic planes?
searcher.add_centroid('CENT3', 'CENT1', 'CENT2')

# searcher.add_distance_constraint('CENT1', 'CENT2', ('==', 7.0), 
#                                     vdw_corrected=False)
# searcher.add_point_plane_distance_constraint('DIST_C', 'CENT1', 'PLANE2', (6.5, 7.5))
# make sure the angle is between -5 and +5 deg?
searcher.add_plane_angle_constraint('ANGLE_C', 'PLANE1', 'PLANE2', (-5, 5))

# Measure the constraints for each substructure found?
searcher.add_plane_angle_measurement('ANGLE', 'PLANE1', 'PLANE2')
searcher.add_distance_constraint('DIST', 'CENT1', 'CENT2', (6.5, 7.5), vdw_corrected=False, type='any')

# How to ensure nothing in-between?

# should I add atom alignment constraints? There's a lot of symmetry in this substructure, and I 
# hope that it won't affect the results, by picking specific carbon atoms to measure angles from..

# find hits
start_time = time.time()
# just get unique hits, compute number of hits per structure in the for loop below.
hits = searcher.search(max_hit_structures=10, max_hits_per_structure=1)
end_time = time.time()

elapsed = end_time - start_time
print('{0:d} Hits found in {1:.2f} seconds.'.format(len(hits), elapsed))
f = open(join(outdir, datafile), 'w')

cwriter = csv.writer(f)
cwriter.writerow(['CSD_NAME', 'UNIT_VOL_A^3', 'CRYSTAL_MOLAR_DENS_MMOL_CM^3', 'SUBSTRUCT_COUNT', 'SUBSTRUCT_DENS_MMOL_CM^3', 
                  'PLANAR_ANGLE', 'PLANAR_ANGLE_STDEV', 'PLANAR_DIST', 'PLANAR_DIST_STDEV'])

# mine just the mofs?
mof_csd = join(csd_subsetdir, 'MOF_subset.gcd')
# want to find the number of hits per structure. (to find density)
for e in ccdc.io.EntryReader(mof_csd):
    print(e.identifier)


for h in hits:
    print(h.identifier)
    local_hits = searcher.search(h.crystal, max_hits_per_structure=10000)
    a, d = zip(*[(s.measurements['ANGLE'], s.constraints['DIST']) for s in local_hits])

    AVG_ANGLE = np.absolute(np.array(a)).mean()
    STD_ANGLE = np.absolute(np.array(a)).std()
    AVG_DIST = np.array(d).mean()
    STD_DIST = np.array(d).std()

    nhits = len(local_hits)
    vol = cm3 / h.crystal.cell_volume * 1000 # how many unit cells fit into a cubic centimetre
    molar_dens = nhits * vol / avo # in mol/cm^3
    cwriter.writerow([h.identifier, h.crystal.cell_volume, vol/avo, nhits, molar_dens, AVG_ANGLE, STD_ANGLE, AVG_DIST, STD_DIST])
    # print('number of hits in unit cell: {0:d}'.format(nhits))
    # print('molar density (mmol/cm^3)  : {0:.3f}'.format(molar_dens))
    # print('mmols per cm^3:              {0:.3f}'.format(vol/avo))
    # print('ANGLE: {0:.2f}'.format(AVG_ANGLE))
    # print('DIST:  {0:.2f}'.format(AVG_DIST))
    writer = ccdc.io.CrystalWriter(join(outdir, '{0:s}.cif'.format(h.identifier)), append=False)
    writer.write(h.crystal)

f.close()

