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
"""Hirshfeld atoms-in-molecules partitioning module."""


import numpy as np

from denspart.stockholder import Stockholder
from denspart.utils import doc_inherit


class Hirshfeld(Stockholder):
    """Hirshfeld Partitioning Scheme."""

    name = "h"

    def __init__(self, coordinates, numbers, pseudo_numbers, dens, grid, proatomdb,
                 procharges=None, min_dens=1.e-100,):
        """
        """
        super(Hirshfeld, self).__init__(coordinates, numbers, pseudo_numbers, dens, grid)

        self._proatomdb = proatomdb
        self._min_dens = min_dens

        # assign initial pro-atom charges
        if procharges is None:
            self._procharges = np.zeros(self.natoms)
        elif len(procharges) == self.natoms:
            self._procharges = np.array(procharges)
        else:
            raise ValueError("Argument pro-charges should have {} length!".format(self.natoms))

    @property
    def procharges(self):
        """Charge of pro-atoms."""
        return self._procharges

    @property
    def proatomdb(self):
        """An instance of `ProAtomDB`."""
        return self._proatomdb

    @property
    def min_dens(self):
        """Minimum value of density."""
        return self._min_dens

    @doc_inherit(Stockholder)
    def compute_proatom_density(self, index):
        # get pro-atom basis expansion
        expansion = self.get_proatom_basis_expansion(index)
        dens = np.zeros(self.grid.size)
        temp = np.zeros(self.grid.size)
        # make linear combination of basis
        for proatom, coeff in expansion.iteritems():
            temp[:] = 0.
            proatom.evaluate(self.coordinates[index], self.grid.points, temp)
            dens += coeff * temp
        assert np.all(dens > 0.)
        return dens

    def get_proatom_basis_expansion(self, index):
        """
        """
        number, charge = self.numbers[index], self.procharges[index]
        # case of integer pro-atom charge
        if abs(charge - int(charge)) < 1.e-6:
            basis = self.proatomdb.get_records(number=number, charge=int(charge))
            if len(basis) == 0:
                raise ValueError("ProAtom with number={} & charge={} dose not "
                                 "exist".format(number, charge))
            elif len(basis) > 1:
                raise ValueError("More than one ProAtom with number={} & charge={} exits, "
                                 "cannot decide what to do!".format(number, charge))
            expansion = {basis[0]: 1.0}
        # case of non-integer pro-atom charge
        else:
            # lower and upper bound charge, coefficient & basis
            lower_q, upper_q = np.floor(charge), np.ceil(charge)
            lower_c, upper_c = upper_q - charge, charge - lower_q
            lower_b = self.proatomdb.get_records(number=number, charge=lower_q)
            upper_b = self.proatomdb.get_records(number=number, charge=upper_q)
            if len(lower_b) == 1 and len(upper_b) == 1:
                expansion = {lower_b[0]: lower_c, upper_b[0]: upper_c}
            elif number - upper_q == 0.:
                expansion = {lower_b[0]: lower_c}
            elif number - upper_q < 0.:
                expansion = {lower_b[0]: 0.}
            else:
                raise ValueError("Cannot make {}-th pro-atom!".format(index))
        return expansion
