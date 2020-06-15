"""This module defines an ASE interface to ORCA 3.0.3
by Ragnar Bjornsson
Based on NWchem interface but simplified.
Only supports energies and gradients (no dipole moments, orbital energies etc.) for now.
For more ORCA-keyword flexibility, method/xc/basis etc. keywords are not used.
Instead two keywords, orcasimpleinput and orcablock are used to define
the ORCA simple-inputline and the ORCA-block input.
This allows for more flexible use of any ORCA method or keyword available in ORCA
instead of hardcoding stuff.
"""
import os

import numpy as np
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError
from ase.io.orca import write_orca
from ase.units import Hartree, Bohr


class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []


class ORCA(FileIOCalculator):
    implemented_properties = ['energy', 'forces']
    command = 'orca PREFIX.inp > PREFIX.out'

    default_parameters = dict(
        charge=0, mult=1,
        task='gradient',
        orcasimpleinput='PBE def2-SVP',
        orcablocks='%scf maxiter 200 end',
    )

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='orca', atoms=None, **kwargs):
        """Construct ORCA-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        # customizing the orca command to use
        if 'orca_command' in kwargs:
            self.command = f'{str(kwargs.get("orca_command"))} PREFIX.inp > PREFIX.out'

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.write(self.label + '.ase')
        with open(self.label + '.inp', 'w') as f:
            f.write(f"! {p.orcasimpleinput} \n")
            f.write(f"{p.orcablocks} \n")
            write_orca(f, atoms, p.charge, p.mult)

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        with open(self.label + '.inp') as f:
            for line in f:
                if line.startswith('geometry'):
                    break
            symbols = []
            positions = []
            for line in f:
                if line.startswith('end'):
                    break
                words = line.split()
                symbols.append(words[0])
                positions.append([float(word) for word in words[1:]])

        self.parameters = Parameters.read(self.label + '.ase')
        self.read_results()

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()

    def read_energy(self):
        """Read Energy from ORCA output file."""
        with open(self.label + '.out', 'r') as f:
            text = f.read()
        lines = iter(text.split('\n'))
        # Energy:
        estring = 'FINAL SINGLE POINT ENERGY'
        for line in lines:
            if estring in line:
                energy = float(line.split()[-1])
                break
        self.results['energy'] = energy * Hartree

    def read_forces(self):
        """Read Forces from ORCA output file."""
        with open(f'{self.label}.engrad', 'r') as file:
            lines = file.readlines()
        getgrad = "no"
        gradients = []
        tempgrad = []
        for i, line in enumerate(lines):
            if line.find('# The current gradient') >= 0:
                getgrad = "yes"
                gradients = []
                tempgrad = []
                continue
            if getgrad == "yes" and "#" not in line:
                grad = line.split()[-1]
                tempgrad.append(float(grad))
                if len(tempgrad) == 3:
                    gradients.append(tempgrad)
                    tempgrad = []
            if '# The at' in line:
                getgrad = "no"
        self.results['forces'] = -np.array(gradients) * Hartree / Bohr
