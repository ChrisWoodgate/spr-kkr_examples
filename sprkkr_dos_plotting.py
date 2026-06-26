'''
sprkkr_dos_plotting.py

Routine(s) for handling SPR-KKR DoS output for pretty plotting in
Python, e.g. using matplotlib.

C. D. Woodgate, University of Bristol, 2026
'''

import numpy as np

ry_to_ev = 13.6057039763

def sprkkr_dos_read(filename, l_max=2, spin_channels=2):
    ''' 
    Routine to read DOS calculation output from SPR-KKR.
        
    NOTE: NOT TESTED FOR ALL CLASSES OF SYSTEM YET
          USE WITH CAUTION
    '''
    
    # Open the DoS file
    dos_file = filename
    
    # Read the data line by line.
    # Typically, SPR-KKR writes the DoS in a weird format with no
    # delimiter, so we will have to break up the strings...
    lines = open(dos_file, 'r')
    data = lines.readlines()

    # Scrape the header of the DoS file for the useful information
    for i, line in enumerate(data):
        if (line[:10].strip() == 'TITLE'):
           title = line[10:].strip()
        elif (line[:10].strip() == 'NT_eff'):
           n_species = int(line[10:].strip())
        elif (line[:10].strip() == 'NE'):
           n_energies = int(line[10:].strip())
        elif (line[:10].strip() == 'EFERMI'):
           e_fermi = float(line[10:].strip())
        elif (line[:10].strip() == 'IT'):
           elements = []
           concentrations = []
           nats = []
           for j in range(i+1, i+n_species+1):
               row = data[j]
               name = row[6:10].strip()
               elements.append(name)
               concentration = float(row[23:30])
               concentrations.append(concentration)
               nat = int(row[34:35])
               nats.append(nat)
           # Note where the actual DoS data will start
           dos_start = j+2
           break

    # Arrays for real components of energies, imaginary components, dos
    # and lm- and species-resolved components
    re_energies = np.zeros(n_energies)
    im_energies = np.zeros(n_energies)
    dos = np.zeros((n_energies))
    lm_components = np.zeros((n_species, n_energies, l_max+1, 
                              spin_channels))
    species_resolved = np.zeros((n_species, n_energies))

    # Loop over energies
    for i in range(n_energies):
        # Loop over species
        all_lines = ''
        for s in range(n_species):
            line = data[dos_start+i*(n_species)+s]
            all_lines += (line.strip())

        re_energies[i] = float(all_lines[0:10])
        im_energies[i] = float(all_lines[10:20])

        start = 20

        for s in range(n_species):
            for j in range(l_max+1):
                for k in range(spin_channels):
                    component = all_lines[start+ j*10 + k*l_max*10 + s*spin_channels*l_max*10:start+ j*10 + k*l_max*10 + s*spin_channels*l_max*10+10]
                    lm_components[s,i,j,k] = float(component)
                    species_resolved[s,i] += float(component)
                    dos[i] += float(component)*nats[s]
        
#            # Handle case of first species (line comes with energy at
#            # start
#            if (s==0):
#                re_energies[i] = float(line[0:10])
#                start = 20
#            # Handle lines below first species
#            else:
#                start = 10
#            # Loop over l values and spins to sum species and total DoS
#            for j in range(l_max):
#                for k in range(spin_channels):
#                    component = line[start+ j*10 + k*l_max*10:start+j*10 + k*l_max*10+10]
#                    if (component == ''):
#                        continue
#                    lm_components[s,i,j,k] = float(component)
#                    species_resolved[s,i] += float(component)
#                    dos[i] += float(component)
#    
    # Convert to units of eV, eV^-1 as relevant
    re_energies = (re_energies-e_fermi)*ry_to_ev
    dos = dos/ry_to_ev
    species_resolved = species_resolved/ry_to_ev

    return re_energies, dos, species_resolved, title, elements, concentrations, n_species 
