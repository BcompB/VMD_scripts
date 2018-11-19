# RMSF to B-factor column in PDB file
#
# Author: Robert Arbon https://github.com/RobertArbon
#
# Python script which takes two input files (MD trajectory and PDB file)
# 	and outputs a PDB file with calculated RMSF values in the final 'B-factor' column.
# In this way, users can select 'color by B-factor' in VMD, and visualise the
# 	selected topology file coloured by RMSF calculated over the MD simulation.

import numpy as np
import mdtraj as md
import pandas as pd
import sys

traj_path = sys.argv[1]
top_path = sys.argv[2]

traj = md.load(traj_path, top=top_path)
top = traj[0]

indices = top.topology.select('protein')

top = top.atom_slice(indices)
print(top)
traj = traj.atom_slice(indices)
print(traj)

traj = traj.superpose(traj)

df, _ = top.topology.to_dataframe()

avg_xyz = np.mean(traj.xyz[:, : , :], axis=0)
msf = 3*np.mean((traj.xyz[:, :, :] - avg_xyz)**2, axis=(0,2))
df['msf'] = msf

by_res = pd.DataFrame(df.groupby(['resSeq'])['msf'].mean())
by_res = by_res.reset_index()
by_res['rmsf'] = np.sqrt(by_res['msf'])

df = pd.merge(left=df, right=by_res, on=['resSeq'])

bfacs = 10*df['rmsf']

print('Mean RMSF  {:4.2f}'.format(np.mean(bfacs)))
print('Min RMSF  {:4.2f}'.format(np.min(bfacs)))
print('Max RMSF  {:4.2f}'.format(np.max(bfacs)))

top.save('rmsf.pdb', bfactors=bfacs)

