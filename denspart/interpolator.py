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
"""Interpolation Module."""


import numpy as np

from horton import CubicSpline, Cell, eval_spline_grid


class AtomInterpolator(object):
    """Atomic Density Interpolator."""

    def __init__(self, rtransform, dens, deriv=None, min_dens=1.e-100):
        """
        """
        self.min_dens = min_dens
        # compute log density
        log_dens = np.zeros_like(dens)
        np.clip(dens, self.min_dens, np.inf, out=log_dens)
        np.log(log_dens, out=log_dens)
        # compute derivative of log density
        if deriv is not None:
            log_deriv = deriv / np.ma.masked_less_equal(dens, 1.e-10)
            log_deriv = np.ma.filled(log_deriv, fill_value=0.0)
        else:
            log_deriv = None
        self.spline = CubicSpline(log_dens, log_deriv, rtransform)

    def evaluate(self, center, points, output=None):
        """
        """
        if output is None:
            output = points.zeros()
        # compute dens on the given points
        eval_spline_grid(self.spline, center, output, points, Cell(None))
        output[output == 0.0] = np.log(self.min_dens)
        np.exp(output, out=output)
        if not (output > 0.0).all():
            raise ValueError("Evaluating a spline on a grid should result in positive values.")
        return output
