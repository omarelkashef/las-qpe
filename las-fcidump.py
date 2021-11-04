import numpy as np
# PySCF imports
from pyscf import gto, scf, lib, mcscf
from pyscf.tools import fcidump
# mrh imports
from mrh.my_pyscf.mcscf.lasscf_o0 import LASSCF

# Define molecule: (H2)_2
xyz = '''H 0.0 0.0 0.0
         H 1.0 0.0 0.0
         H 0.2 3.9 0.1
         H 1.159166 4.1 -0.1'''
mol = gto.M (atom = xyz, basis = 'sto-3g', output='h4_sto3g.log',
    verbose=lib.logger.DEBUG)

# Do RHF
mf = scf.RHF(mol).run()

# Create LASSCF object
las = LASSCF(mf, (2,2), (2,2), spin_sub=(1,1))

# Localize the chosen fragment active spaces
frag_atom_list = ((0,1), (2,3))
mf.mo_coeff = las.localize_init_guess(frag_atom_list, mf.mo_coeff)

print ("Ncore: ",las.ncore, "Ncas: ",las.ncas, "Ncas_sub", las.ncas_sub)
for i, sub in enumerate(las.ncas_sub):
    if i == 0:
        fcidump.from_mo(mol, 'fcidump_las_h4_sub{}'.format(i),mf.mo_coeff[:,:las.ncore+sub])
    else:
        fcidump.from_mo(mol, 'fcidump_las_h4_sub{}'.format(i), np.hstack((mf.mo_coeff[:,:las.ncore], mf.mo_coeff[:,las.ncas_sub[i-1]:las.ncas_sub[i-1]+sub])))

