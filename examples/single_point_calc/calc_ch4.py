"""
Example:
Calculating single point energy and forces for methane
"""

from ase.build import molecule
from ase.calculators.orca import ORCA

# change for your local path, unless python can see the executable simply like this
_orca_command = "orca"

if __name__ == '__main__':
    calc = ORCA(label="orca",
                orca_command=_orca_command,
                charge=0,
                mult=1,
                task='gradient',
                orcasimpleinput='engrad RHF revPBE def2-TZVP def2/J D3BJ slowconv kdiis',
                orcablocks=
                f"%scf Convergence tight \n maxiter 500 \n end"
                )

    atoms = molecule("CH4")
    atoms.set_calculator(calc)

    print("calculating RKS energy and forces for methane structure from ase.build.molecule:")
    energy = atoms.get_potential_energy()
    print(f"energy: {energy:10.6f} eV")
    print("forces:", atoms.get_forces())
