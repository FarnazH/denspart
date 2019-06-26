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


import os
import numpy as np

from glob import glob

from horton import IOData, RadialGrid, ExpRTransform, BeckeMolGrid
from horton import ProAtomDB as HORTONProAtomDB
from denspart.hirshfeld import Hirshfeld, HirshfeldI
from denspart.proatomdb import ProAtomDB, ProAtomRecord
from denspart.interpolator import AtomInterpolator


def test_hirshfeld_charges_h2o():
    # get pro-atom database
    fpath = os.path.join(os.path.dirname(__file__),
                         "cached/atom_???_???_hf_sto3g.fchk")
    proatomdb = HORTONProAtomDB.from_files(glob(fpath))

    # convert to new ProAtomDB
    records = []
    for number in proatomdb.get_numbers():
        for charge in proatomdb.get_charges(number):
            record = proatomdb.get_record(number, charge)
            rtransform = proatomdb.get_rgrid(number).rtransform
            inter = AtomInterpolator(rtransform, record.rho, record.deriv)
            records.append(ProAtomRecord(inter, number=number, charge=charge))
    proatomdb = ProAtomDB(records)

    # get molecule
    fpath = os.path.join(os.path.dirname(__file__),
                         "cached/water_sto3g_hf_g03.fchk")
    mol = IOData.from_file(fpath)
    dm_full = mol.get_dm_full()

    # get molecular grid
    rtf = ExpRTransform(5e-4, 2e1, 120)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        (rgrid, 110), random_rotate=False, mode="discard")
    # compute molecular density
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid.points)

    # check computed charges
    part = Hirshfeld(mol.coordinates, mol.numbers, mol.pseudo_numbers, moldens, grid, proatomdb)
    part.run()
    expected = np.array([-0.246171541212, 0.123092011074, 0.123079530138])
    assert abs(part.charges - expected).max() < 1.e-4


def test_hirshfeldi_charges_h2o():
    # get pro-atom database
    fpath = os.path.join(os.path.dirname(__file__),
                         "cached/atom_???_???_hf_sto3g.fchk")
    proatomdb = HORTONProAtomDB.from_files(glob(fpath))

    # convert to new ProAtomDB
    records = []
    for number in proatomdb.get_numbers():
        for charge in proatomdb.get_charges(number):
            record = proatomdb.get_record(number, charge)
            rtransform = proatomdb.get_rgrid(number).rtransform
            inter = AtomInterpolator(rtransform, record.rho, record.deriv)
            records.append(ProAtomRecord(inter, number=number, charge=charge))
    proatomdb = ProAtomDB(records)

    # get molecule
    fpath = os.path.join(os.path.dirname(__file__),
                         "cached/water_sto3g_hf_g03.fchk")
    mol = IOData.from_file(fpath)
    dm_full = mol.get_dm_full()

    # get molecular grid
    rtf = ExpRTransform(5e-4, 2e1, 120)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        (rgrid, 110), random_rotate=False, mode="discard")
    # compute molecular density
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid.points)

    # check computed charges
    part = HirshfeldI(mol.coordinates, mol.numbers, mol.pseudo_numbers, moldens, grid, proatomdb)
    part.run()
    expected = np.array([-0.4214, 0.2107, 0.2107])
    print(part.charges)
    print(expected)
    print(abs(part.charges - expected).max())
    print(abs(part.charges - expected).max() < 1.e-3)
    assert abs(part.charges - expected).max() < 1.e-3
