"""
This module contains a chamber prmtop class that will read in all
parameters and allow users to manipulate that data and write a new
prmtop object.

Copyright (C) 2010 - 2015  Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
   
You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""
from __future__ import division

from chemistry.amber._amberparm import AmberParm
from chemistry.constants import NTYPES, NATYP, IFBOX, TINY, NATOM
from chemistry.exceptions import AmberParmError
from chemistry.topologyobjects import (UreyBradley, Improper, Cmap, BondType,
                                       ImproperType, CmapType)
from math import sqrt

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ChamberParm(AmberParm):
    """
    Chamber Topology (parm7 format) class. Gives low, and some high, level
    access to topology data or interact with some of the high-level classes
    comprising the system topology and parameters.  The ChamberParm class uses
    the same attributes that the AmberParm class uses, and only the ones unique
    to ChamberParm will be shown below.

    Parameters
    ----------
    prm_name : str=None
        If provided, this file is parsed and the data structures will be loaded
        from the data in this file
    rst7_name : str=None
        If provided, the coordinates and unit cell dimensions from the provided
        Amber inpcrd/restart file will be loaded into the molecule

    Attributes
    ----------
    In addition to the attributes listed for `AmberParm`, the following
    attributes are available:

    LJ_14_radius : list(float)
        The same as LJ_radius, except specific for 1-4 nonbonded parameters,
        which may differ in the CHARMM force field
    LJ_14_depth : list(float)
        The same as LJ_depth, except specific for 1-4 nonbonded parameters,
        which may differ in the CHARMM force field
    urey_bradleys : TrackedList(UreyBradley)
        List of Urey-Bradley terms between two atoms in a valence angle
    impropers : TrackedList(Improper)
        List of CHARMM-style improper torsions
    cmaps : TrackedList(Cmap)
        List of coupled-torsion correction map parameters
    urey_bradley_types : TrackedList(UreyBradleyType)
        List of parameters defining the Urey-Bradley terms
    improper_types : TrackedList(Improper)
        List of parameters defining the Improper terms
    cmap_types : TrackedList(CmapType)
        List of parameters defining the CMAP terms
    chamber : bool=True
        On ChamberParm instances, this is always True to indicate that it is a
        CHAMBER-style topology file
    amoeba : bool=False
        On ChamberParm instances, this is always False to indicate that it is
        not an AMOEBA-style topology file
    has_cmap : bool
        True if CMAP parameters are present in this system; False otherwise
    """

    solvent_residues = ('WAT', 'TIP3', 'HOH', 'TIP4', 'TIP5', 'SPCE', 'SPC')

    #===================================================

    def initialize_topology(self, rst7_name=None):
        """
        Initializes topology data structures, like the list of atoms, bonds,
        etc., after the topology file has been read. The following methods are
        called:
        """
        self.LJ_14_radius = []
        self.LJ_14_depth = []
        AmberParm.initialize_topology(self, rst7_name)
      
    #===================================================

    def load_pointers(self):
        """
        Loads the data in POINTERS section into a pointers dictionary with each
        key being the pointer name according to http://ambermd.org/formats.html
        """
        AmberParm.load_pointers(self)
        # Other pointers
        nub, nubtypes = self.parm_data['CHARMM_UREY_BRADLEY_COUNT'][:2]
        self.pointers['NUB'] = nub
        self.pointers['NUBTYPES'] = nubtypes
        self.pointers['NIMPHI'] = self.parm_data['CHARMM_NUM_IMPROPERS'][0]
        self.pointers['NIMPRTYPES'] = self.parm_data['CHARMM_NUM_IMPR_TYPES'][0]
        # If CMAP is not present, don't load the pointers
        if self.has_cmap:
            self.pointers['CMAP'] = self.parm_data['CHARMM_CMAP_COUNT'][0]
            self.pointers['CMAP_TYPES'] = self.parm_data['CHARMM_CMAP_COUNT'][1]

    #===================================================

    def load_structure(self):
        """ 
        Loads all of the topology instance variables. This is necessary if we
        actually want to modify the topological layout of our system
        (like deleting atoms)
        """
        super(ChamberParm, self).load_structure()
        self._load_urey_brad_info()
        self._load_improper_info()
        self._load_cmap_info()
        super(ChamberParm, self).unchange()

    #===================================================

    @classmethod
    def load_from_structure(cls, struct):
        """
        Take a Structure instance and initialize a ChamberParm instance from
        that data.

        Parameters
        ----------
        struct : Structure
            The input structure from which to construct a ChamberParm instance
        """
        inst = struct.copy(cls, split_dihedrals=True)
        inst.pointers = {}
        inst.LJ_types = {}
        nbfixes = inst.atoms.assign_nbidx_from_types()
        # Fill the Lennard-Jones arrays/dicts
        ntyp = 0
        for atom in inst.atoms:
            inst.LJ_types[atom.type] = atom.nb_idx
            ntyp = max(ntyp, atom.nb_idx)
        inst.LJ_radius = [0 for i in xrange(ntyp)]
        inst.LJ_depth = [0 for i in xrange(ntyp)]
        inst.LJ_14_radius = [0 for i in xrange(ntyp)]
        inst.LJ_14_depth = [0 for i in xrange(ntyp)]
        for atom in inst.atoms:
            inst.LJ_radius[atom.nb_idx-1] = atom.atom_type.rmin
            inst.LJ_depth[atom.nb_idx-1] = atom.atom_type.epsilon
            inst.LJ_14_radius[atom.nb_idx-1] = atom.atom_type.rmin_14
            inst.LJ_14_depth[atom.nb_idx-1] = atom.atom_type.epsilon_14
        inst._add_standard_flags()
        inst.pointers['NATOM'] = len(inst.atoms)
        inst.parm_data['POINTERS'][NATOM] = len(inst.atoms)
        if struct.box is None:
            inst.parm_data['POINTERS'][IFBOX] = 0
            inst.pointers['IFBOX'] = 0
        elif (abs(struct.box[3] - 90) > TINY or abs(struct.box[4] - 90) > TINY
                or abs(struct.box[5] - 90) > TINY):
            inst.parm_data['POINTERS'][IFBOX] = 2
            inst.pointers['IFBOX'] = 2
            inst.parm_data['BOX_DIMENSIONS'] = [struct.box[3]] + struct.box[:3]
        else:
            inst.parm_data['POINTERS'][IFBOX] = 1
            inst.pointers['IFBOX'] = 1
            inst.parm_data['BOX_DIMENSIONS'] = [90] + struct.box[:3]
        try:
            coords = []
            for atom in struct.atoms:
                coords.extend([atom.xx, atom.xy, atom.xz])
        except AttributeError:
            raise
        else:
            inst.load_coordinates(coords)
        inst.remake_parm()
        inst._set_nonbonded_tables(nbfixes)

        return inst

    #===================================================

    def remake_parm(self):
        """
        Fills :attr:`parm_data` from the data in the parameter and topology
        arrays (e.g., :attr:`atoms`, :attr:`bonds`, :attr:`bond_types`, ...)
        """
        # Get rid of terms containing deleted atoms and empty residues
        self.prune_empty_terms()
        self.residues.prune()
        self.rediscover_molecules()

        # Transfer information from the topology lists 
        self._xfer_atom_info()
        self._xfer_residue_info()
        self._xfer_bond_info()
        self._xfer_angle_info()
        self._xfer_dihedral_info()
        self._xfer_urey_bradley_properties()
        self._xfer_improper_properties()
        self._xfer_cmap_properties()
        # Load the pointers dict
        self.load_pointers()
        # Mark atom list as unchanged
        super(ChamberParm, self).unchange()

    #===================================================

    def fill_LJ(self):
        """
        Fills the LJ_radius, LJ_depth arrays and LJ_types dictionary with data
        from LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF sections of the prmtop
        files, by undoing the canonical combining rules.
        """
        AmberParm.fill_LJ(self)

        acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        ntypes = self.pointers['NTYPES']

        self.LJ_14_radius = []  # empty LJ_radii so it can be re-filled
        self.LJ_14_depth = []   # empty LJ_depths so it can be re-filled
        one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

        for i in xrange(ntypes):
            lj_index = self.parm_data["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if acoef[lj_index] < 1.0e-6:
                self.LJ_14_radius.append(0)
                self.LJ_14_depth.append(0)
            else:
                factor = 2 * acoef[lj_index] / bcoef[lj_index]
                self.LJ_14_radius.append(pow(factor, one_sixth) * 0.5)
                self.LJ_14_depth.append(bcoef[lj_index] / 2 / factor)

    #===================================================

    def recalculate_LJ(self):
        """
        Takes the values of the LJ_radius and LJ_depth arrays and recalculates
        the LENNARD_JONES_A/BCOEF topology sections from the canonical
        combining rules.
        """
        AmberParm.recalculate_LJ(self)
        ntypes = self.pointers['NTYPES']
        acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        for i in xrange(ntypes):
            for j in xrange(i, ntypes):
                index = self.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = self.LJ_14_radius[i] + self.LJ_14_radius[j]
                wdij = sqrt(self.LJ_14_depth[i] * self.LJ_14_depth[j])
                acoef[index] = wdij * rij ** 12
                bcoef[index] = 2 * wdij * rij**6

    #===================================================

    @property
    def chamber(self):
        return True
   
    @property
    def amoeba(self):
        return False

    @property
    def has_cmap(self):
        return len(self.cmaps) > 0 or 'CHARMM_CMAP_COUNT' in self.parm_data

    #===========  PRIVATE INSTANCE METHODS  ============

    def _load_urey_brad_info(self):
        """ Loads the Urey-Bradley types and array """
        del self.urey_bradleys[:]
        del self.urey_bradley_types[:]
        for k, req in zip(self.parm_data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT'],
                          self.parm_data['CHARMM_UREY_BRADLEY_EQUIL_VALUE']):
            self.urey_bradley_types.append(
                    BondType(k, req, self.urey_bradley_types)
            )
        ulist = self.parm_data['CHARMM_UREY_BRADLEY']
        for i in xrange(0, 3*self.pointers['NUB'], 3):
            self.urey_bradleys.append(
                    UreyBradley(self.atoms[ulist[i  ]-1],
                                self.atoms[ulist[i+1]-1],
                                self.urey_bradley_types[ulist[i+2]-1])
            )

    #===================================================

    def _load_improper_info(self):
        """ Loads the CHARMM Improper types and array """
        del self.impropers[:]
        del self.improper_types[:]
        for k, eq in zip(self.parm_data['CHARMM_IMPROPER_FORCE_CONSTANT'],
                         self.parm_data['CHARMM_IMPROPER_PHASE']):
            self.improper_types.append(
                    ImproperType(k, eq, self.improper_types)
            )
        ilist = self.parm_data['CHARMM_IMPROPERS']
        for i in xrange(0, 5*self.pointers['NIMPHI'], 5):
            self.impropers.append(
                    Improper(self.atoms[ilist[i  ]-1],
                             self.atoms[ilist[i+1]-1],
                             self.atoms[ilist[i+2]-1],
                             self.atoms[ilist[i+3]-1],
                             self.improper_types[ilist[i+4]-1])
            )

    #===================================================

    def _load_cmap_info(self):
        """ Loads the CHARMM CMAP types and array """
        if not self.has_cmap: return
        del self.cmaps[:]
        del self.cmap_types[:]
        for i in xrange(self.pointers['CMAP_TYPES']):
            resolution = self.parm_data['CHARMM_CMAP_RESOLUTION'][i]
            grid = self.parm_data['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            cmts = self.parm_comments['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            self.cmap_types.append(
                    CmapType(resolution, grid, cmts, list=self.cmap_types)
            )
        clist = self.parm_data['CHARMM_CMAP_INDEX']
        for i in xrange(0, 6*self.pointers['CMAP'], 6):
            self.cmaps.append(
                    Cmap(self.atoms[clist[i  ]-1], self.atoms[clist[i+1]-1],
                         self.atoms[clist[i+2]-1], self.atoms[clist[i+3]-1],
                         self.atoms[clist[i+4]-1],
                         self.cmap_types[clist[i+5]-1])
            )

    #===================================================

    def _check_section_lengths(self):
        """ Make sure each section has the necessary number of entries """
        super(ChamberParm, self)._check_section_lengths()

        def check_length(key, length, required=True):
            if not required and not key in self.parm_data: return
            if len(self.parm_data[key]) != length:
                raise AmberParmError('FLAG %s has %d elements; expected %d' %
                                     (key, len(self.parm_data[key]), length))

        check_length('CHARMM_UREY_BRADLEY_COUNT', 2)
        check_length('CHARMM_UREY_BRADLEY', self.pointers['NUB']*3)
        check_length('CHARMM_UREY_BRADLEY_FORCE_CONSTANT',
                     self.pointers['NUBTYPES'])
        check_length('CHARMM_UREY_BRADLEY_EQUIL_VALUE',
                     self.pointers['NUBTYPES'])
        check_length('CHARMM_NUM_IMPROPERS', 1)
        check_length('CHARMM_IMPROPERS', self.pointers['NIMPHI']*5)
        check_length('CHARMM_NUM_IMPR_TYPES', 1)
        check_length('CHARMM_IMPROPER_FORCE_CONSTANT',
                     self.pointers['NIMPRTYPES'])
        check_length('CHARMM_IMPROPER_PHASE', self.pointers['NIMPRTYPES'])

        ntypes = self.pointers['NTYPES']
        check_length('LENNARD_JONES_14_ACOEF', ntypes*(ntypes+1)//2)
        check_length('LENNARD_JONES_14_BCOEF', ntypes*(ntypes+1)//2)
        if self.has_cmap:
            check_length('CHARMM_CMAP_COUNT', 2)
            check_length('CHARMM_CMAP_RESOLUTION', self.pointers['CMAP_TYPES'])
            for i in range(self.pointers['CMAP_TYPES']):
                res = self.parm_data['CHARMM_CMAP_RESOLUTION'][i]
                check_length('CHARMM_CMAP_PARAMETER_%02d' % (i+1), res*res)

    #===================================================

    def _xfer_urey_bradley_properties(self):
        """
        Sets the various topology file section data from the Urey-Bradley arrays
        """
        data = self.parm_data
        for urey_type in self.urey_bradley_types:
            urey_type.used = False
        for urey in self.urey_bradleys:
            urey.type.used = True
        self.urey_bradley_types.prune_unused()
        data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT'] = \
                [type.k for type in self.urey_bradley_types]
        data['CHARMM_UREY_BRADLEY_EQUIL_VALUE'] = \
                [type.req for type in self.urey_bradley_types]
        data['CHARMM_UREY_BRADLEY'] = bond_array = []
        for urey in self.urey_bradleys:
            bond_array.extend([urey.atom1.idx+1, urey.atom2.idx+1,
                               urey.type.idx+1])
        nub = len(self.urey_bradleys)
        nubt = len(self.urey_bradley_types)
        data['CHARMM_UREY_BRADLEY_COUNT'] = [nub, nubt]
        self.pointers['NUB'] = nub
        self.pointers['NUBTYPES'] = nubt

    #===================================================

    def _xfer_improper_properties(self):
        """ Sets the topology file section data from the improper arrays """
        data = self.parm_data
        for improper_type in self.improper_types:
            improper_type.used = False
        for improper in self.impropers:
            improper.type.used = True
        self.improper_types.prune_unused()
        data['CHARMM_IMPROPER_FORCE_CONSTANT'] = \
                [type.psi_k for type in self.improper_types]
        data['CHARMM_IMPROPER_PHASE'] = \
                [type.psi_eq for type in self.improper_types]
        data['CHARMM_IMPROPERS'] = improper_array = []
        for imp in self.impropers:
            improper_array.extend([imp.atom1.idx+1, imp.atom2.idx+1,
                                   imp.atom3.idx+1, imp.atom4.idx+1,
                                   imp.type.idx+1])
        data['CHARMM_NUM_IMPROPERS'] = [len(self.impropers)]
        data['CHARMM_NUM_IMPR_TYPES'] = [len(self.improper_types)]
        self.pointers['NIMPHI'] = len(self.impropers)
        self.pointers['NIMPRTYPES'] = len(self.improper_types)

    #===================================================

    def _xfer_cmap_properties(self):
        """ Sets the topology file section data from the cmap arrays """
        if not self.has_cmap: return
        data = self.parm_data
        for ct in self.cmap_types:
            ct.used = False
        for cmap in self.cmaps:
            cmap.type.used = True
        self.cmap_types.prune_unused()
        # If we have no cmaps, delete all remnants of the CMAP terms in the
        # prmtop and bail out
        if len(self.cmaps) == 0:
            # We have deleted all cmaps. Get rid of them from the parm file and
            # bail out. This is probably pretty unlikely, though...
            self.delete_flag('CHARMM_CMAP_COUNT')
            self.delete_flag('CHARMM_CMAP_RESOLUTION')
            for flag in self.flag_list:
                if flag.startswith('CHARMM_CMAP_PARAMETER'):
                    self.delete_flag(flag)
            del self.pointers['CMAP']
            del self.pointers['CMAP_TYPES']
            return
        # All of our CMAP types are in different topology file sections. We need
        # to delete all of the CHARMM_CMAP_PARAMETER_XX sections and then
        # recreate them with the correct size and comments.  The comments have
        # been stored in the CMAP types themselves to prevent them from being
        # lost. We will also assume that the Fortran format we're going to use
        # is the same for all CMAP types, so just pull it from
        # CHARMM_CMAP_PARAMETER_01 (or fall back to 8(F9.5))
        try:
            fmt = str(self.formats['CHARMM_CMAP_PARAMETER_01'])
        except KeyError:
            fmt = '8(F9.5)'
        flags_to_delete = []
        for flag in self.flag_list:
            if flag.startswith('CHARMM_CMAP_PARAMETER'):
                flags_to_delete.append(flag)
        for flag in flags_to_delete:
            self.delete_flag(flag)
        # Now add them back
        for i, ct in enumerate(self.cmap_types):
            self.add_flag('CHARMM_CMAP_PARAMETER_%02d' % (i+1), fmt,
                          data=ct.grid, comments=ct.comments)
        # Now do the CMAP_INDEX section
        data['CHARMM_CMAP_INDEX'] = cmap_array = []
        for cm in self.cmaps:
            cmap_array.extend([cm.atom1.idx+1, cm.atom2.idx+1,
                               cm.atom3.idx+1, cm.atom4.idx+1,
                               cm.atom5.idx+1, cm.type.idx+1])
        data['CHARMM_CMAP_COUNT'] = [len(self.cmaps), len(self.cmap_types)]
        data['CHARMM_CMAP_RESOLUTION'] = \
                    [ct.resolution for ct in self.cmap_types]
        self.pointers['CMAP'] = len(self.cmaps)
        self.pointers['CMAP_TYPES'] = len(self.cmap_types)

    #===================================================

    def _add_standard_flags(self):
        """ Adds all of the standard flags to the parm_data array """
        self.set_version()
        self.add_flag('CTITLE', '20a4', num_items=0)
        self.add_flag('POINTERS', '10I8', num_items=31)
        self.add_flag('FORCE_FIELD_TYPE', 'i2,a78', num_items=0)
        self.add_flag('ATOM_NAME', '20a4', num_items=0)
        self.add_flag('CHARGE', '3E24.16', num_items=0,
                comments=['Atomic charge multiplied by sqrt(332.0716D0) '
                         '(CCELEC)'])
        self.add_flag('ATOMIC_NUMBER', '10I8', num_items=0)
        self.add_flag('MASS', '5E16.8', num_items=0)
        self.add_flag('ATOM_TYPE_INDEX', '10I8', num_items=0)
        self.add_flag('NUMBER_EXCLUDED_ATOMS', '10I8', num_items=0)
        self.add_flag('NONBONDED_PARM_INDEX', '10I8', num_items=0)
        self.add_flag('RESIDUE_LABEL', '20a4', num_items=0)
        self.add_flag('RESIDUE_POINTER', '10I8', num_items=0)
        self.add_flag('BOND_FORCE_CONSTANT', '5E16.8', num_items=0)
        self.add_flag('BOND_EQUIL_VALUE', '5E16.8', num_items=0)
        self.add_flag('ANGLE_FORCE_CONSTANT', '5E16.8', num_items=0)
        self.add_flag('ANGLE_EQUIL_VALUE', '3E25.17', num_items=0)
        self.add_flag('CHARMM_UREY_BRADLEY_COUNT', '2I8', num_items=2,
                comments=['V(ub) = K_ub(r_ik - R_ub)**2',
                          'Number of Urey Bradley terms and types'])
        self.add_flag('CHARMM_UREY_BRADLEY', '10I8', num_items=0,
                comments=['List of the two atoms and its parameter index',
                          'in each UB term: i,k,index'])
        self.add_flag('CHARMM_UREY_BRADLEY_FORCE_CONSTANT', '5E16.8',
                num_items=0, comments=['K_ub: kcal/mol/A**2'])
        self.add_flag('CHARMM_UREY_BRADLEY_EQUIL_VALUE', '5E16.8', num_items=0,
                comments=['r_ub: A'])
        self.add_flag('DIHEDRAL_FORCE_CONSTANT', '5E16.8', num_items=0)
        self.add_flag('DIHEDRAL_PERIODICITY', '5E16.8', num_items=0)
        self.add_flag('DIHEDRAL_PHASE', '5E16.8', num_items=0)
        self.add_flag('SCEE_SCALE_FACTOR', '5E16.8', num_items=0)
        self.add_flag('SCNB_SCALE_FACTOR', '5E16.8', num_items=0)
        self.add_flag('CHARMM_NUM_IMPROPERS', '10I8', num_items=0,
                comments=['Number of terms contributing to the',
                          'quadratic four atom improper energy term:',
                          'V(improper) = K_psi(psi - psi_0)**2'])
        self.add_flag('CHARMM_IMPROPERS', '10I8', num_items=0,
                comments=['List of the four atoms in each improper term',
                          'i,j,k,l,index  i,j,k,l,index',
                          'where index is into the following two lists:',
                          'CHARMM_IMPROPER_{FORCE_CONSTANT,IMPROPER_PHASE}'])
        self.add_flag('CHARMM_NUM_IMPR_TYPES', '1I8', num_items=1,
                comments=['Number of unique parameters contributing to the',
                          'quadratic four atom improper energy term'])
        self.add_flag('CHARMM_IMPROPER_FORCE_CONSTANT', '5E16.8', num_items=0,
                comments=['K_psi: kcal/mole/rad**2'])
        self.add_flag('CHARMM_IMPROPER_PHASE', '5E16.8', num_items=0,
                comments=['psi: degrees'])
        natyp = self.pointers['NATYP'] = self.parm_data['POINTERS'][NATYP] = 1
        self.add_flag('SOLTY', '5E16.8', num_items=natyp)
        self.add_flag('LENNARD_JONES_ACOEF', '3E24.16', num_items=0)
        self.add_flag('LENNARD_JONES_BCOEF', '3E24.16', num_items=0)
        self.add_flag('LENNARD_JONES_14_ACOEF', '3E24.16', num_items=0)
        self.add_flag('LENNARD_JONES_14_BCOEF', '3E24.16', num_items=0)
        self.add_flag('BONDS_INC_HYDROGEN', '10I8', num_items=0)
        self.add_flag('BONDS_WITHOUT_HYDROGEN', '10I8', num_items=0)
        self.add_flag('ANGLES_INC_HYDROGEN', '10I8', num_items=0)
        self.add_flag('ANGLES_WITHOUT_HYDROGEN', '10I8', num_items=0)
        self.add_flag('DIHEDRALS_INC_HYDROGEN', '10I8', num_items=0)
        self.add_flag('DIHEDRALS_WITHOUT_HYDROGEN', '10I8', num_items=0)
        self.add_flag('EXCLUDED_ATOMS_LIST', '10I8', num_items=0)
        self.add_flag('HBOND_ACOEF', '5E16.8', num_items=0)
        self.add_flag('HBOND_BCOEF', '5E16.8', num_items=0)
        self.add_flag('HBCUT', '5E16.8', num_items=0)
        self.add_flag('AMBER_ATOM_TYPE', '20a4', num_items=0)
        self.add_flag('TREE_CHAIN_CLASSIFICATION', '20a4', num_items=0)
        self.add_flag('JOIN_ARRAY', '10I8', num_items=0)
        self.add_flag('IROTAT', '10I8', num_items=0)
        if self.has_cmap:
            self.add_flag('CHARMM_CMAP_COUNT', '2I8', num_items=2,
                    comments=['Number of CMAP terms, number of unique CMAP '
                              'parameters'])
            self.add_flag('CHARMM_CMAP_RESOLUTION', '20I4', num_items=0,
                    comments=['Number of steps along each phi/psi CMAP axis',
                              'for each CMAP_PARAMETER grid'])
            self.add_flag('CHARMM_CMAP_INDEX', '6I8', num_items=0,
                    comments=['Atom index i,j,k,l,m of the cross term',
                              'and then pointer to CHARMM_CMAP_PARAMETER_n'])
        if self.box is not None:
            self.add_flag('SOLVENT_POINTERS', '3I8', num_items=3)
            self.add_flag('ATOMS_PER_MOLECULE', '10I8', num_items=0)
            self.add_flag('BOX_DIMENSIONS', '5E16.8', num_items=4)
        self.add_flag('RADIUS_SET', '1a80', num_items=1)
        self.add_flag('RADII', '5E16.8', num_items=0)
        self.add_flag('SCREEN', '5E16.8', num_items=0)
        self.add_flag('IPOL', '1I8', num_items=1)

    #===================================================

    def _set_nonbonded_tables(self, nbfixes=None):
        """
        Sets the tables of Lennard-Jones nonbonded interaction pairs
        """
        data = self.parm_data
        ntypes = data['POINTERS'][NTYPES]
        ntypes2 = ntypes * ntypes
        # Set up the index lookup tables (not a unique solution)
        data['NONBONDED_PARM_INDEX'] = [0 for i in xrange(ntypes2)]
        holder = [0 for i in xrange(ntypes2)]
        idx = 0
        for i in xrange(ntypes):
            for j in xrange(i+1):
                idx += 1
                holder[ntypes*i+j] = holder[ntypes*j+i] = idx
        idx = 0
        for i in xrange(ntypes):
            for j in xrange(ntypes):
                data['NONBONDED_PARM_INDEX'][idx] = \
                            holder[ntypes*i+j]
                idx += 1
        nttyp = ntypes * (ntypes + 1) // 2
        # Now build the Lennard-Jones arrays
        data['LENNARD_JONES_14_ACOEF'] = [0 for i in xrange(nttyp)]
        data['LENNARD_JONES_14_BCOEF'] = [0 for i in xrange(nttyp)]
        data['LENNARD_JONES_ACOEF'] = [0 for i in xrange(nttyp)]
        data['LENNARD_JONES_BCOEF'] = [0 for i in xrange(nttyp)]
        self.recalculate_LJ()
        # Now make any NBFIX modifications we had
        if nbfixes is not None:
            for i, fix in enumerate(nbfixes):
                for terms in fix:
                    j, rmin, eps, rmin14, eps14 = terms
                    i, j = min(i, j-1), max(i, j-1)
                    eps = abs(eps)
                    eps14 = abs(eps14)
                    idx = data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                    data['LENNARD_JONES_ACOEF'][idx] = eps * rmin**12
                    data['LENNARD_JONES_BCOEF'][idx] = 2 * eps * rmin**6
                    data['LENNARD_JONES_14_ACOEF'][idx] = eps14 * rmin14**12
                    data['LENNARD_JONES_14_BCOEF'][idx] = 2 * eps14 * rmin14**6

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

charmm_solvent = ('WAT', 'TIP3', 'HOH', 'TIP4', 'TIP5', 'SPCE', 'SPC')

def ConvertFromPSF(struct, params, title=''):
    """
    This function instantiates a ChamberParm instance from a data structure
    instantiated by a CHARMM PSF.

    Parameters
    ----------
    struct : Structure
        The structure object (typically loaded from a PSF file)
    params : CharmmParameterSet
        The parameter set describing the parameters of the input struct
    title : str=''
        The title to assign to the topology file

    Returns
    -------
    ChamberParm
        ChamberParm instance with all parameters loaded
    """
    # Make sure all atom types are strings, not integers
    int_starting = str(struct.atoms[0].atom_type) != struct.atoms[0].type
    if int_starting:
        for atom in struct.atoms:
            atom.type = str(atom.atom_type)
    parm = ChamberParm.load_from_structure(struct)
    parm.parm_data['FORCE_FIELD_TYPE'] = fftype = []
    for pset in params.parametersets:
        if 'CHARMM' not in pset: # needed to trigger "charmm_active"...
            pset = 'CHARMM: %s' % pset
        fftype.extend([len(params.parametersets), pset])
    if not fftype:
        fftype.extend([1, 'CHARMM force field: No FF information parsed...'])

    # Convert atom types back to integers if that's how they started
    if int_starting:
        for atom in struct.atoms:
            atom.type = int(atom.atom_type)
    return parm
