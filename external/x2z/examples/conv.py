import automol

with open('ccc.xyz', 'r') as f:
    xyz_str = f.read()

xyz = automol.geom.from_string(xyz_str)
for x in xyz:
    print(x)
print(xyz)
