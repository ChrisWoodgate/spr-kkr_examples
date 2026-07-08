# SPR-KKR calculations on fcc Cu
#
# C. D. Woodgate, 2026

def main():

    # Import the relevant modules
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from ase import Atoms
    from ase.build import bulk

    # I was using a local version of ase2sprkkr
    # You should not need the below two lines if you have pip installed the package
    #import sys
    #sys.path.append('/path/to/ase2sprkkr/src')

    from ase2sprkkr.sprkkr.calculator import SPRKKR
    from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms
    from ase.eos import EquationOfState
    from ase.units import Bohr,Rydberg,kJ,kB,fs,Hartree,mol,kcal

    #=================================================#
    # Step 1: Self-consistent field (SCF) calculation #
    #=================================================#

    # lattice parameter
    a_lat = 3.62

    # Make the directory we need (check if it exists first).
    # (This means the SCF doesn't clutter the main directory)
    directory = f'./scf'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Call 'bulk' to build me a bulk fcc crystal
    atom = bulk('Cu', 'fcc', a=a_lat)

    # Now call ase2sprkkr to let me work on the atoms object
    # in SPR-KKR mode. This allows specification of things like
    # partial lattice site occupancies (for the CPA).
    # (Don't actually need to do it here, but included for completeness.)
    atoms=SPRKKRAtoms.promote_ase_atoms(atom)

    # Now we set up the ASE calculator
    # Specify how you run SPR-KKR here, along with controlling your PATH, etc.
    calculator = SPRKKR(atoms=atoms,mpi=['mpirun','-np','96'], directory=directory)

    # And the necessary parameters
    calculator.input_parameters.set(NKTAB=500) # Control k-point sampling. 250 is default.
    calculator.input_parameters.set(VXC='VWN') # VWN (LDA) XC functional
    calculator.input_parameters.set(MODE='SREL') # Scalar relativistic mode
    calculator.input_parameters.set(NL=3) # l_max=2 means 3 l-channels
    calculator.input_parameters.set(NE=32) # 32-point semicircular energy contour
    calculator.input_parameters.ENERGY.GRID=5 # Semicircular energy contour
    calculator.input_parameters.SCF.MIX=0.05 # Gentle mixing paramter to ensure convergence
    calculator.input_parameters.SCF.TOL=1e-6 # Tight tolerance on SCF cycle
    calculator.input_parameters.SCF.NITER=200 # Allow for lots of iterations if needed
    calculator.input_parameters.ENERGY.ImE=0.005 # Closest approach to real axis
    
    # Run the calculator
    out=calculator.calculate()

    #========================================================#
    # Step 2: Electronic density of states (DOS) calculation #
    #========================================================#

    Ef = out.last_iteration['energy']['EF'].to_dict()

    top = Ef + 0.5
    bottom = Ef - 1.0

    converged_pot = out.potential_filename.split('/')[-1]

    calculator = out.calculator

    # Make the directory we need (check if it exists first).
    directory = f'./dos/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    os.system(f'cp ./scf/' + converged_pot + f' ./dos/')

    calculator.input_parameters='DOS'
    calculator.input_parameters.set(NKTAB=5000) # Control k-point sampling. 250 is default.
    calculator.input_parameters.ENERGY.EMAX=top
    calculator.input_parameters.ENERGY.EMIN=bottom
    calculator.input_parameters.ENERGY.ImE=0.01 # Closest approach to real axis
    calculator.input_parameters.set(NE=1501)
    calculator.input_parameters.ENERGY.GRID=3 # Parallel energy contour

    calculator.directory = directory

    dos = calculator.calculate(potential=directory+converged_pot, print_output='info')

    #===================================================#
    # Step 3: Bloch spectral function (BSF) calculation #
    #===================================================#

    calculator = out.calculator

    # Make the directory we need (check if it exists first).
    directory = f'./blochsf/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    os.system(f'cp ./scf/' + converged_pot + f' ./blochsf/')

    calculator.input_parameters='bsfek'
    calculator.input_parameters.TASK.KPATH=1 # Parallel energy contour
    calculator.input_parameters.TASK.NK=1001 # Parallel energy contour
    calculator.input_parameters.set(NKTAB=5000) # Control k-point sampling. 250 is default.
    calculator.input_parameters.ENERGY.EMAX=top
    calculator.input_parameters.ENERGY.EMIN=bottom
    calculator.input_parameters.ENERGY.ImE=0.001 # Closest approach to real axis
    calculator.input_parameters.set(NE=1501)
    calculator.input_parameters.ENERGY.GRID=3 # Parallel energy contour

    calculator.directory = directory

    bsf = calculator.calculate(potential=directory+converged_pot, print_output='info')

# If we run this script, execute 'main'
if __name__ == "__main__":
    main()
