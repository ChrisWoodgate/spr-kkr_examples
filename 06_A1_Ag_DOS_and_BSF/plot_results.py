import sys
sys.path.append('..')
import sprkkr_dos_plotting
import sprkkr_bsf_plotting
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

# Font size etc.
plt.rc('font',**{'size'   : 16,
             'weight' :'normal'})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Parse the BSF data
bsfdata = sprkkr_bsf_plotting.parse_bsfdata('blochsf/_Ag_BSFEK_BLOCHSF.bsf')

# Parse the DOS data
re_energies, dos, species_resolved, title, elements, concentrations, n_species = sprkkr_dos_plotting.sprkkr_dos_read('dos/_Ag_DOS_DOS.dos', l_max=2, spin_channels=2)

# Make a figure and grab my axes
fig,axes = plt.subplots(1,2, sharey=True,gridspec_kw={'width_ratios': [2.4, 1]})
fig.set_size_inches(7,4)
ax1 = axes[0]
ax2 = axes[1]

# Maximum BSF value
Imax = np.max(np.abs(np.sum(bsfdata.I,axis=0)))

# Plot the bsf
# You may need to play around with vmax to make this look nice
ax1.pcolormesh(bsfdata.k,(bsfdata.E-bsfdata.EFERMI)*sprkkr_bsf_plotting.Ry,np.sum(bsfdata.I,axis=0),shading='gouraud',cmap='Greys',vmin=0.0,vmax=Imax/3, rasterized=True)

# Top and bottom energies for limits
Emin = np.min((bsfdata.E-bsfdata.EFERMI)*sprkkr_bsf_plotting.Ry)
Emax = np.max((bsfdata.E-bsfdata.EFERMI)*sprkkr_bsf_plotting.Ry)
ax1.set_ylim([-12.5,5])    

# Add the ticks, lines and special point labels
ax1.set_xticks(np.insert(bsfdata.k[0,[X-1 for X in bsfdata.INDKDIR]],0,0))

# Uncomment the line below if using KPATH=1
#ax1.set_xticklabels([r'$X$', r'$\Gamma$', r'$L$', r'$W$', r'$K$', r'$\Gamma/L$', r'$U$', r'$X$', r'$W$', r'$U$'])
# Uncomment the line below if using KPATH=2
ax1.set_xticklabels([r'$X$', r'$\Gamma$', r'$L$', r'$W$', r'$K$', r'$\Gamma$'])

for index in bsfdata.INDKDIR:
    ax1.plot([bsfdata.k[0,index-1],bsfdata.k[0,index-1]],[-12.5,Emax],color='black',lw=0.5)

# Add the Fermi level
ax1.plot([0,np.max(bsfdata.k)],[0,0],color='black',lw=0.5, linestyle='dashed')

# Some labels
ax1.set_title(r'A1 Ag')
ax1.set_ylabel(r'$E-E_{\rm F}$ (eV)')

# Plot the DOS
ax2.hlines(0.0, 0.0, 12.0 ,color='black',lw=0.5, linestyle='dashed')
ax2.set_xlabel(r'DOS (eV$^{-1}$)')
ax2.set_xlim(0.0, 4.0)
for i in range(0, n_species,2):
    ax2.plot(species_resolved[i,:], re_energies, alpha=1.0, label=elements[i], color='tab:gray')
ax2.legend(fontsize=11, loc='lower right')

plt.tight_layout()
plt.savefig('./Ag.pdf', dpi=300)

nelec = simpson(dos[0:1002], x=re_energies[0:1002])
print('Sanity check: Found ', nelec, ' electrons up to Fermi energy by integrating raw DOS data')
