# -*- coding: utf-8 -*-
# DensPart: Density-Based Atoms-in-Molecules Partitioning Package.
#
# Copyright (C) 2019 The DensPart Development Team
#
# This file is part of DensPart.
#
# DensPart is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# DensPart is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
"""Stockholder atoms-in-molecules partitioning module."""


import numpy as np


class Stockholder(object):
    """Stockholder Partitioning Scheme."""

    def __init__(self, coordinates, numbers, pseudo_numbers, dens, grid):
        """Initialize class.

        Parameters
        ----------
        coordinates: np.ndarray(M, 3)
            Cartesian coordinates of M atoms in molecule.
        numbers : np.ndarray(M,)
            Atomic number of M atoms in molecule.
        pseudo_numbers : np.ndarray(M,)
            Pseudo atomic number of M atoms in molecule.
        dens : np.ndarray(N,)
            Molecular density evaluated on N grid points.
        grid : instance of Grid
            Instance of molecular grid class.

        """
        # check size of arguments
        if grid.size != len(dens):
            raise ValueError("Arguments dens & grid represent different number of points.")
        if len(numbers) != len(coordinates):
            raise ValueError("Arguments numbers & coordinates represent different number of atoms.")
        # assign attributes
        self._coordinates = coordinates
        self._numbers = numbers
        self._pseudo_numbers = pseudo_numbers
        self._dens = dens
        self._grid = grid
        # assign empty arrays
        self._charges = np.zeros(self.natoms)
        self._weights = np.zeros((self.natoms, self.grid.size))
        self._populations = np.zeros(self.natoms)
        self._prodens = np.zeros(self.grid.size)
        self._procharges = np.zeros(self.natoms)

    @property
    def coordinates(self):
        """Cartesian coordinates of atoms in molecule."""
        return self._coordinates

    @property
    def numbers(self):
        """Atomic number of atoms in molecule."""
        return self._numbers

    @property
    def natoms(self):
        """Number of atoms in the molecule."""
        return len(self._numbers)

    @property
    def pseudo_numbers(self):
        """Pseudo atomic number of atoms in molecule."""
        return self._pseudo_numbers

    @property
    def grid(self):
        """Instance of grid class."""
        return self._grid

    @property
    def dens(self):
        """Molecular density evaluated on the grid points."""
        return self._dens

    @property
    def prodens(self):
        """Pro-molecular density evaluated on the grid points."""
        return self._prodens

    @property
    def weights(self):
        """Atomic weights evaluated on the grid points."""
        return self._weights

    @property
    def charges(self):
        """Charge of atoms in molecule."""
        return self._charges

    @property
    def populations(self):
        """Population of atoms in molecule."""
        return self._populations

    def run(self):
        """Run the atoms-in-molecule partitioning."""
        # compute reference pro-atom densities
        for index in range(self.natoms):
            self._weights[index, :] = self.compute_proatom_density(index)
        # compute stockholder atoms
        self.compute_stockholder_partition()

    def compute_proatom_density(self, index):
        """Compute density of pro-atom with the given index.

        Parameters
        ----------
        index : int
            Index of pro-atom.

        Returns
        -------
        out : np.ndarray(N,)
            Density of pro-atom evaluated on N grid points.

        """
        raise NotImplementedError

    def compute_stockholder_partition(self):
        """Compute stockholder atomic weights, charges and populations."""
        # compute pro-molecule density & atomic weights
        self._prodens = np.sum(self._weights, axis=0)
        self._prodens[self._prodens < self.min_dens] = self.min_dens
        self._weights /= self._prodens
        # compute atomic population & charge
        self._populations = np.apply_along_axis(self.grid.integrate, 1, self.weights * self.dens)
        self._populations += self.numbers - self.pseudo_numbers
        self._charges = self._numbers - self._populations
