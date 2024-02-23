import pymatgen.core as pmg


li2o = pmg.Structure.from_spacegroup(
    "Fm-3m", pmg.Lattice.cubic(2), ["Li", "O"], [[0.25, 0.25, 0.25], [0, 0, 0]]
)

#print(li2o)
print(li2o.get_all_neighbors(10))

#print(type(li2o[2]))
