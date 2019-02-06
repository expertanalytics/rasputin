from rasputin import triangulate_dem

"""
Most minimal example
"""
vs = triangulate_dem.PointVector([(1,0,0), (0,1,0), (0,0,0), (0.25,0.25,1)])

points, faces = triangulate_dem.lindstrom_turk_by_ratio(vs, 2.0)

print("\n Points:")
for p in points:
    print('  ', p)

print("\n Faces:")
for f in faces:
    print('  ', f)

normals = triangulate_dem.surface_normals(points, faces)

print("\n Normals:")
for n in normals:
    print('  ', n)
