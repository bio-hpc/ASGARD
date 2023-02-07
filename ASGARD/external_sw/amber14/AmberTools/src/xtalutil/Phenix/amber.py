#! /usr/bin/env phenix.python
from libtbx import group_args
import sys, os
import iotbx.pdb
import argparse
from scitbx.array_family import flex
import phenix_amber_interface as ext
import scitbx.restraints
#~ import amber_adaptbx_ext as PAI
import time


master_phil_str = """
  use_amber = False
    .type = bool
  topology_file_name = None
    .type = path
  coordinate_file_name = None
    .type = path
"""

class geometry_manager(object):

  def __init__(self,
        sites_cart=None,
        energy_components=flex.double([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),
        gradients=None,
        number_of_restraints=0,
        gradients_factory=flex.double,
        mdgx_structs=None):
    self.sites_cart = sites_cart
    self.energy_components = energy_components
    self.gradients_factory = gradients_factory
    self.number_of_restraints=number_of_restraints
    self.mdgx_structs=mdgx_structs

  def energies_sites(self, compute_gradients=False):
    #Convert flex arrays to C arrays
    sites_cart_c=ext.ExtractVec(self.sites_cart.as_double())
    gradients_c=ext.ExtractVec(flex.double(self.sites_cart.size() * 3, 0))
    energy_components_c=ext.ExtractVec(self.energy_components)

    # Call c++ interface to call mdgx to calculate new gradients and target
    ext.callMdgx(sites_cart_c, gradients_c, energy_components_c, self.mdgx_structs)
    if (compute_gradients) :
      gradients = self.gradients_factory(gradients_c) * -1
    else :
      gradients = self.gradients_factory(
        flex.double(self.sites_cart.size() * 3,0))
    result = energies(
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=None,
      gradients_factory=None,
      normalization=False)
    result.number_of_restraints = self.number_of_restraints
    result.residual_sum = float(energy_components_c[0])
    result.energy_components = list(energy_components_c)
    result.finalize_target_and_gradients()
    #result.show()
    return result

class energies (scitbx.restraints.energies) :
  def __init__ (self, *args, **kwds) :
    scitbx.restraints.energies.__init__(self, *args, **kwds)
    self.energy_components = None

  def show(self):
    print "\n\n"
    print "Amber_total_energy: %7.6f"     %(self.residual_sum)
    print "  bonds (n=%d): %7.6f"       %(self.energy_components[6],
      self.energy_components[1])
    print "  angles (n=%d): %7.6f"      %(self.energy_components[7],
      self.energy_components[2])
    print "  dihedrals (n=%d): %7.6f"     %(self.energy_components[8],
      self.energy_components[3])
    print "  electrostatics: %7.6f"     %(self.energy_components[4])
    print "  vanderWaals: %7.6f"      %(self.energy_components[5])
    print "\n\n"
    return 0

  def get_rmsd_gradient(self):
    return self.energy_components[9]

def print_sites_cart(sites_cart):
        for atom in sites_cart:
                print("%8.3f%8.3f%8.3f"%(atom[0], atom[1], atom[2]))

def get_amber_structs (prmtop, ambcrd):
        return ext.uform(prmtop, ambcrd)

def run(pdb,prmtop, crd):

  #===================================================================#
  #                                                                   #
  #  BEFORE C++                                                       #
  #                                                                   #
  #===================================================================#

  #file i/o
  pdb_file = os.path.abspath(pdb)
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  pdb_atoms = pdb_inp.atoms_with_labels()
  symm = pdb_inp.crystal_symmetry()
  xray_structure = pdb_inp.xray_structure_simple(enable_scattering_type_unknown=True)


  #     initiate flex arrays for coordinates, gradients, energy
  sites_cart=xray_structure.sites_cart()
  gradients=flex.double(len(sites_cart)*3)
  target=flex.double([6.7,1.0,2.0,3.0,4.0,5.0,0.0,0.0,0.0,0.0])
  print "Number of atom sites: %d " %sites_cart.size()
  print "\nGradients and target BEFORE C call:"
  print list(gradients[1:10])
  print target[0]

  #===================================================================#
  #                                                                   #
  #  CALL C++                                                         #
  #                                                                   #
  #===================================================================#

  U=ext.uform(prmtop, crd)

  #Convert flex arrays to C arrays
  sites_cart_c=ext.ExtractVec(sites_cart.as_double())
  gradients_c=ext.ExtractVec(gradients)
  target_c=ext.ExtractVec(target)

  # Call c++ interface to call mdgx to calculate new gradients and target
  ext.callMdgx(sites_cart_c, gradients_c, target_c, U)

  # Convert back into python types (eg. into flex arrays for phenix to use)
  gradients=flex.vec3_double(gradients_c)*-1

  target= flex.double(target_c)

  #===================================================================#
  #                                                                   #
  #  AFTER C++                                                        #
  #                                                                   #
  #===================================================================#


  print "\nGradients and target AFTER C call:"
  print list(gradients[0:10])
  print target[0]
  print target[9]

  print "\n"
  print "Amber_total_energy: %7.6f"             %(target[0])
  print "  bonds (n= ): %7.6f"                  %(target[1])
  print "  angles (n= ): %7.6f"                         %(target[2])
  print "  dihedrals (n= ): %7.6f"              %(target[3])
  print "  electrostatics: %7.6f"               %(target[4])
  print "  vanderWaals: %7.6f"                  %(target[5])
  print "\n\n"

  return 0

if __name__ == "__main__" :
        parser = argparse.ArgumentParser()
        parser.add_argument("pdb", help="name of pdb file")
        parser.add_argument("prmtop", help="name of topology file")
        parser.add_argument("crd", help="name of coordinate file")
        args = parser.parse_args()
        run(args.pdb,args.prmtop, args.crd)
