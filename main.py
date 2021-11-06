#!/usr/bin/env python
'''
https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/substructure_searching.html

Substructure searching with geometric measurements

'''
import time
import ccdc.search
import ccdc.io


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
searcher.add_point_plane_distance_constraint('DIST_C', 'CENT1', 'PLANE2', (6.5, 7.5))
# make sure the angle is between -5 and +5 deg?
searcher.add_plane_angle_constraint('ANGLE_C', 'PLANE1', 'PLANE2', (-5, 5))

# Measure the constraints for each substructure found?
searcher.add_plane_angle_measurement('ANGLE', 'PLANE1', 'PLANE2')
searcher.add_distance_measurement('DIST', 'CENT1', 'CENT2')

# should I add atom alignment constraints? There's a lot of symmetry in this substructure, and I hope that it won't affect the results, by picking specific carbon atoms to measure angles from..


# find hits

hits = searcher.search()

for h in hits:
    print('ANGLE: {0:.2f}'.format(h.measurements['ANGLE']))
    print('DIST:  {0:.2f}'.format(h.measurements['DIST']))


#sub_id = substructure_search.add_substructure(testosterone_substructure)
#hits = substructure_search.search()
#print (len(hits))
#print(len(substructure_search.search(max_hit_structures=4)))
# to search an individual structure:
# substructure_search.search(structure)

