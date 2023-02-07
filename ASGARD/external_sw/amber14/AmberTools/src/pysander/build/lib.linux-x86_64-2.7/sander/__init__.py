from chemistry.amber import AmberParm, Rst7
from chemistry import unit as u
import tempfile
import sys as _sys
try:
    import numpy as _np
except ImportError:
    _np = None

__all__ = ['InputOptions', 'QmInputOptions', 'setup', 'cleanup', 'pme_input',
           'gas_input', 'natom', 'energy_forces', 'set_positions', 'set_box',
           'is_setup', 'EnergyTerms']

from array import array as _array
try:
    if _sys.version_info[0] > 2:
        exec('from . import pysander as _pys')
    else:
        exec('import pysander as _pys')
except ImportError:
    raise ImportError('Could not import the compiled Python-sander interface. '
                      'Make sure you have the Python development libraries '
                      'installed and try rebuilding the serial installation '
                      'of AMBER.')

# Accept unit-ized input from OpenMM if that is the source of the input data
try:
    from simtk.unit import is_quantity, angstroms
except ImportError:
    is_quantity = lambda *args, **kwargs: False

# If set to True, units are applied to the resulting output. Otherwise,
# everything is left unitless (in AKMA units)
APPLY_UNITS = False

# Add some of the pysander members directly to the sander namespace
InputOptions = _pys.InputOptions
QmInputOptions = _pys.QmInputOptions
EnergyTerms = _pys.EnergyTerms
cleanup = _pys.cleanup
pme_input = _pys.pme_input
gas_input = _pys.gas_input
natom = _pys.natom
is_setup = _pys.is_setup

# To help with dimensional analysis handling
def _strip_units(obj):
    """
    Strips units from the object and returns its value in the AKMA unit system.
    If it is a scalar, the original object is returned unchanged
    """
    if u.is_quantity(obj):
        return obj.value_in_unit_system(u.akma_unit_system)

def _strip_units_from_struct(struct):
    """
    Strips all units from all members of the struct-like class
    """
    for attr in dir(struct):
        val = _strip_units(getattr(struct, attr))
        setattr(struct, attr, val)
    return struct

def _apply_units_to_struct(struct, unit):
    """ Applies the given unit to all members of the struct """
    for attr in dir(struct):
        val = getattr(struct, attr)
        setattr(struct, attr, val*unit)
    return struct

# For Python3 compatibility
try:
    basestring
except NameError:
    basestring = str

# Ideally we'd use contextlib.contextmanager and turn the function into a
# context manager directly as a generator. However, contextlib was added in
# Python 2.5, and we desire to keep Python 2.4 support while still allowing the
# use of the context manager. As a result, we turn setup into a class (yuck)
# with __enter__ and __exit__ methods that implement the context manager for
# Python 2.5+, or do nothing in Python 2.4. It behaves the same way that the
# function would.

class setup(object):
    """ Sets up a sander calculation. This supports acting like a context
    manager such that the cleanup routine is called upon exiting the context. It
    can also be called as a standalone function.

    Parameters
    ----------
    prmtop : AmberParm or str
        Name of the prmtop file to use to set up the calculation or an AmberParm
        instance
    coordinates : list/iterable or str
        list of coordinates or name of the inpcrd file
    box : list/iterable or None
        list of 3 box lengths and 3 box angles. Can be None if no box is
        required or if the box will be read in from the coordinates
    mm_options : InputOptions
        struct with sander options
    qm_options : QmInputOptions (optional)
        struct with the QM options in sander QM/MM calculations

    Examples
    --------

    The following are equivalent invocations which each make sure that the
    sander data structures are cleaned up afterwards

    >>> with sander.setup("prmtop", inpcrd.coords, inpcrd.box, mm_options):
    ...     e, f = sander.energy_forces()
    ... 
    >>> sander.is_setup()
    False

    Without a context manager (e.g., Python 2.4 and earlier), which ensures that
    cleanup is done

    >>> try:
    ...     sander.setup("prmtop", inpcrd.coords, inpcrd.box, mm_options)
    ...     e, f = sander.energy_forces()
    ... finally:
    ...     if sander.is_setup():
    ...         sander.cleanup()
    ... 
    >>> sander.is_setup()
    False

    If you want the sander system to stay set up when the current function ends
    and persist until another function call, simply do not execute a cleanup in
    a `finally` clause and do not use a context manager

    >>> sander.setup("prmtop", inpcrd.coords, inpcrd.box, mm_options)
    >>> e, f = sander.energy_forces()
    >>> sander.is_setup()
    True
    """

    def __init__(self, prmtop, coordinates, box, mm_options, qm_options=None):
        # Handle the case where the coordinates are actually a restart file
        if isinstance(coordinates, basestring):
            # This is a restart file name. Parse it and make sure the coordinates
            # and box
            rst = Rst7.open(coordinates)
            try:
                coordinates = rst.coordinates.tolist()
            except AttributeError:
                coordinates = rst.coordinates
            if rst.hasbox and not box:
                try:
                    box = rst.box.tolist()
                except AttributeError:
                    box = rst.box

        # Convert from numpy arrays to regular arrays
        if hasattr(coordinates, 'tolist'): # works for numpy.ndarray and array.array
            coordinates = coordinates.tolist()
        if hasattr(box, 'tolist'):
            box = box.tolist()
        if not box:
            box = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        try:
            box = [float(x) for x in box]
        except TypeError:
            raise TypeError('box must be an iterable with 6 numerical elements')
        if len(box) != 6:
            raise ValueError('box must have 6 elements')

        # Check if the prmtop is an AmberParm instance or not. If it is, write out a
        # temporary prmtop file
        if isinstance(prmtop, AmberParm):
            parm = tempfile.mktemp(suffix='.parm7')
            prmtop.write_parm(parm)
        elif not isinstance(prmtop, basestring):
            raise TypeError('prmtop must be an AmberParm or string')
        else:
            parm = prmtop

        # Error checking
        if mm_options.ifqnt != 0 and qm_options is None:
            raise ValueError("qm_options must be provided if QM/MM is requested")

        # Call the setup routine
        if qm_options is None:
            _pys.setup(parm, coordinates, box, mm_options)
        else:
            _pys.setup(parm, coordinates, box, mm_options, qm_options)

    def __enter__(self):
        """ Nothing needs to be done here """
        pass

    def __exit__(self, *args, **kwargs):
        """ Make sure that sander is cleaned up """
        if _pys.is_setup(): _pys.cleanup()

def qm_input():
    """
    Returns a populated set of QM input options. Only here for consistency with
    gas_input and pme_input -- QmInputOptions can be instantiated directly
    """
    return QmInputOptions()

def set_positions(positions):
    """
    Sets the particle positions of the active system from the passed list of
    positions. Supports both lists, numpy.ndarray and numpy.ndarray objects

    Parameters
    ----------
    positions : array of float
        The atomic positions. They can have units of length. They can have the
        shapes (natom*3,) or (natom, 3)
    """
    if is_quantity(positions):
        positions = positions.value_in_unit(angstroms)
    # Common input types will have an natom x 3 shape. I can call "flatten" on
    # numpy arrays to solve this quickly, but in cases where the coordinates
    # are given as a list (or tuple) of Vec3's (or tuples), this requires
    # separate handling
    try:
        positions = positions.flatten()
    except AttributeError:
        natom = _pys.natom()
        if len(positions) == natom:
            p = positions
            positions = [0.0 for i in range(natom*3)]
            for i, x in enumerate(p):
                i3 = i * 3
                try:
                    positions[i3], positions[i3+1], positions[i3+2] = x
                except ValueError:
                    raise ValueError('Expected iterable with shape (natom, 3)')
        else:
            if len(positions) != natom * 3:
                raise ValueError('Positions array must have shape (natom, 3) '
                                 'or (natom)')
    if hasattr(positions, 'tolist'): # works for array.array and numpy.ndarray
        return _pys.set_positions(positions.tolist())
    return _pys.set_positions(positions)

def get_positions(as_numpy=False):
    """ Returns the current atomic positions loaded in the sander API

    Parameters
    ----------
    as_numpy : bool, optional
        If True, the positions will be returned as a natom*3-length numpy array.
        If False (default), it will be returned as a natom*3-length Python list.

    Returns
    -------
    positions : array of float
        The atomic positions as a list (or numpy array if requested). If
        sander.APPLY_UNITS is True, the return object will be a Quantity with
        the units chemistry.unit.angstroms
    """
    global APPLY_UNITS
    positions = _pys.get_positions()
    if as_numpy:
        positions = np.asarray(positions)
    if APPLY_UNITS:
        return u.Quantity(positions, u.angstrom)
    return positions

def energy_forces(as_numpy=False):
    """
    Returns the energies and forces of the current conformation with the current
    Hamiltonian.

    Parameters
    ----------
    as_numpy : bool, optional
        If True, the forces will be returned as a natom*3-length numpy array. If
        False (default), they will be returned as a natom*3-length Python list.

    Returns
    -------
    energy, forces : EnergyTerms, array of float
        The energies returned in an EnergyTerms container, and the forces are
        returned as a natom*3-length list (or numpy array if requested). If
        sander.APPLY_UNITS is True, the energies will have the units
        kilocalories_per_mole applied, and forces will have the units
        kilocalories_per_mole/u.angstroms
    """
    global APPLY_UNITS
    e, f = _pys.energy_forces()
    if as_numpy:
        f = np.asarray(f)
    if APPLY_UNITS:
        return (_apply_units_to_struct(e, u.kilocalories_per_mole),
                u.Quantity(f, u.kilocalories_per_mole/u.angstroms))
    return e, f

def set_box(a, b, c, alpha, beta, gamma):
    """ Sets the unit cell dimensions for the current system

    Parameters
    ----------
    a : float
        Length of the first unit cell vector (can be a unit.Quantity object
        with dimension length). Unitless input is assumed to be in Angstroms
    b : float
        Length of the second unit cell vector (can be a unit.Quantity object
        with dimension length). Unitless input is assumed to be in Angstroms
    c : float
        Length of the third unit cell vector (can be a unit.Quantity object
        with dimension length). Unitless input is assumed to be in Angstroms
    alpha : float
        Angle between vectors b and c (can be a unit.Quantity object with
        dimension angle). Unitless input is assumed to be in Degrees.
    beta : float
        Angle between vectors a and c (can be a unit.Quantity object with
        dimension angle). Unitless input is assumed to be in Degrees.
    gamma : float
        Angle between vectors a and b (can be a unit.Quantity object with
        dimension angle). Unitless input is assumed to be in Degrees.
    """
    if u.is_quantity(a): a = a.value_in_unit(u.angstroms)
    if u.is_quantity(b): b = b.value_in_unit(u.angstroms)
    if u.is_quantity(c): c = c.value_in_unit(u.angstroms)
    if u.is_quantity(alpha): alpha = alpha.value_in_unit(u.degrees)
    if u.is_quantity(beta): beta = beta.value_in_unit(u.degrees)
    if u.is_quantity(gamma): gamma = gamma.value_in_unit(u.degrees)
    _pys.set_box(a, b, c, alpha, beta, gamma)
