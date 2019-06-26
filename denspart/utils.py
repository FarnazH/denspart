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
"""Utility Function Module."""



def doc_inherit(base_class, base_method=None):
    """Inherit method docstring from another class and/or method.

    Parameters
    ----------
    base_class : str
        Name of base class.
    base_method
        Name of method in the base class.

    """

    def decorator(method):
        """Overwrite method docstring."""
        # check whether the method exists
        if base_method is None:
            overridden = getattr(base_class, method.__name__, None)
        else:
            overridden = getattr(base_class, base_method, None)
        if overridden is None:
            raise AttributeError('Can\'t find method \'%s\' in base class.')
        # change docstring
        method.__doc__ = overridden.__doc__
        return method

    return decorator
