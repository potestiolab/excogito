import mdtraj
import os
print("Sample script for converting GROMACS xtc files in xyz removing hydrogen atoms")

xtc_path = input("insert path to XTC file\n")
gro_path = input("insert path to GRO file\n")
xyz_filename = input("insert path to output XYZ file\n")
full_traj = mdtraj.load_xtc(xtc_path.strip(),top=gro_path.strip())
full_traj_topology = full_traj.topology
print("full trajectory has", len(full_traj_topology.select("all")), "atoms")
# getting no_h atoms
no_h = full_traj_topology.select('type != H')
n_heavy_traj = len(no_h)
print("full trajectory has", n_heavy_traj, "heavy atoms") 
# writing down xyz
mdt_tr_heavy = mdtraj.load_xtc(xtc_path, top=gro_path, atom_indices = list(no_h))
print("number of frames to convert to xyz: ", len(mdt_tr_heavy[:].xyz))
print(mdt_tr_heavy[:].xyz)
with mdtraj.formats.XYZTrajectoryFile(xyz_filename , 'w') as f:
    f.write(mdt_tr_heavy[:].xyz*10)
