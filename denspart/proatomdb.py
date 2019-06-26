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
"""Pro-atom Database Module."""


import numpy as np


class ProAtomRecord(object):
    """Pro-atom Record."""

    def __init__(self, interpolator, **tags):
        self.interpolator = interpolator
        self.tags = tags

    def __getitem__(self, attr):
        value = self.tags.get(attr, None)
        return value

    def __eq__(self, other):
        return self.tags == other.tags

    def __ne__(self, other):
        return not self.__eq__(other)

    def evaluate(self, center, points, output=None):
        return self.interpolator.evaluate(center, points, output)


class ProAtomDB(object):
    """Database of Pro-atoms."""

    def __init__(self, records):
        """Initialize the class.

        Parameters
        ----------
        records : sequence of ProAtomRecord instances
            A sequence of ProAtomRecord instances.

        """
        self._records = records

    def __len__(self):
        """Number of pro-atom records in the database."""
        return len(self._records)

    def get_records(self, **tags):
        selected = []
        for record in self._records:
            checks = [np.array_equal(record.tags.get(k), v) for k, v in tags.items()]
            if all(checks):
                selected.append(record)
        return selected
