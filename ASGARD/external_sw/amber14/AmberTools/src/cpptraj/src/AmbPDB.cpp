/*! \file AmbPDB.cpp
    \brief Potential replacement for ambpdb using cpptraj file routines.
    \author Daniel R. Roe
 */
#include <cstdio> // stdin, fileno
#include <unistd.h> // isatty
#include "CpptrajStdio.h"
#include "ParmFile.h"
#include "Trajin_Single.h"
#include "Trajout_Single.h"
#include "StringRoutines.h"
#include "Traj_AmberRestart.h"
#define VERSION_STRING "V15.00"

static void Help(const char* prgname, bool showAdditional) {
  mprinterr("\nUsage: %s -p 'Top' -c 'Coords' [Additional Options]\n"
              "       %s -p 'Top' < 'AmberRestart' [Additional Options]\n"
              "    -p 'Top'       Topology file (default: prmtop).\n"
              "    -c 'Coords'    Coordinate file.\n"
              "    'AmberRestart' Amber restart file from STDIN.\n"
              "  PDB is written to STDOUT.\n", prgname, prgname);
  if (showAdditional) {
    mprinterr(
            "  Additional Options:\n"
            "    -tit <TITLE>  Write a REMARK record containing TITLE.\n"
            "                      (default: use prmtop title)\n"
            "    -aatm         Left-justified Amber atom names.\n"
            "    -bres         Brookhaven Residue names (HIE->HIS, etc.).\n"
            "    -ctr          Center molecule on (0,0,0).\n"
            "    -noter        Do not write TER records.\n"
            "    -ext          Use PRMTOP extended PDB info, if present.\n"
//            "    -ene <FLOAT>  Define H-bond energy cutoff for FIRST.\n"
//            "    -bin          The coordinate file is in binary form.\n"
            "    -offset <INT> Add offset to residue numbers.\n"
            "  Options for alternate output format (give only one of these):\n"
            "    -pqr          PQR (MEAD) format with charges and radii.\n"
//            "    -sas          PQR with 1.4 added to atom radii.\n"
            "    -mol2         TRIPOS MOL2 format.\n");
//            "    -bnd          list bonds from the PRMTOP.\n"
//            "    -atm          Mike Connolly surface/volume format.\n"
//            "    -first        Add REMARKs for input to FIRST.\n"
  } else
    mprinterr("  Use '-h' or '--help' to see additional options.\n");
  mprinterr("\n");
}

static bool Unsupported(std::string const& arg) {
  if (arg == "-ext" || arg == "-ene" || arg == "-bin" ||
      arg == "-sas" || arg == "-bnd" || arg == "-atm" ||
      arg == "-first") return true;
  return false;
}

// ----- M A I N ---------------------------------------------------------------
int main(int argc, char** argv) {
  SetWorldSilent(true); // No STDOUT output from cpptraj routines.
  mprinterr("| ambpdb (C++) Version %s\n", VERSION_STRING);
  std::string topname, crdname, title, aatm(" pdbatom"), bres, pqr;
  TrajectoryFile::TrajFormatType fmt = TrajectoryFile::PDBFILE;
  bool ctr_origin = false;
  bool noTER = false;
  bool useExtendedInfo = false;
  int res_offset = 0;
  int debug = 0;
  int numSoloArgs = 0;
  for (int i = 1; i < argc; ++i) {
    std::string arg( argv[i] );
    if (arg == "-p" && i+1 != argc && topname.empty()) // Topology
      topname = std::string( argv[++i] );
    else if (arg == "-c" && i+1 != argc && crdname.empty()) // Coords
      crdname = std::string( argv[++i] );
    else if (arg == "-tit" && i+1 != argc && title.empty()) // Title
      title = " title " + std::string( argv[++i] );
    else if (arg == "-offset" && i+1 != argc) // Residue # offset
      res_offset = convertToInteger( argv[++i] );
    else if ((arg == "-d" || arg == "--debug") && i+1 != argc) // Debug level
      debug = convertToInteger( argv[++i] );
    else if (arg == "-h" || arg == "--help") {// Help
      Help(argv[0], true);
      return 0;
    } else if (arg == "-aatm") // Amber atom names
      aatm.clear();
    else if (arg == "-bres") // PDB residue names
      bres.assign(" pdbres");
    else if (arg == "-ext") // Use extended PDB info from Topology
      useExtendedInfo = true;
    else if (arg == "-ctr")  // Center on origin
      ctr_origin = true;
    else if (arg == "-noter") // No TER cards
      noTER = true;
    else if (arg == "-pqr") { // Charge/Radii in occ/bfactor cols
      pqr.assign(" dumpq");
      ++numSoloArgs;
    } else if (arg == "-mol2") { // output as mol2
      fmt = TrajectoryFile::MOL2FILE;
      ++numSoloArgs;
    } else if (Unsupported(arg)) {
      mprinterr("Error: Option '%s' is not yet supported.\n\n", arg.c_str());
      return 1;
    } else {
      mprinterr("Error: Unrecognized option '%s'\n", arg.c_str());
      Help(argv[0], false);
      return 1;
    }
  }
  // Check command line for errors.
  if (topname.empty()) topname.assign("prmtop");
  if (crdname.empty())
    mprinterr("| Reading Amber restart from STDIN\n");
  if (numSoloArgs > 1) {
    mprinterr("Error: Only one alternate output format option may be specified (found %i)\n",
              numSoloArgs);
    Help(argv[0], true);
    return 1;
  }
  if (debug > 0) {
    mprinterr("Warning: debug is %i; debug info will be written to STDOUT.\n", debug);
    SetWorldSilent(false);
  }
  // Topology
  ParmFile pfile;
  Topology parm;
  if (pfile.ReadTopology(parm, topname, debug)) {
    if (topname == "prmtop") Help(argv[0], false);
    return 1;
  }
  parm.IncreaseFrames( 1 );
  if (noTER)
    parm.ClearMoleculeInfo();
  if (!useExtendedInfo)
    parm.ResetPDBinfo();
  if (res_offset != 0)
    for (int r = 0; r < parm.Nres(); r++)
      parm.SetRes(r).SetOriginalNum( parm.Res(r).OriginalResNum() + res_offset );
  ArgList trajArgs;
  // Input coords
  Frame TrajFrame;
  // TODO: Set coord info in parm?
  if (!crdname.empty()) {
    Trajin_Single trajin;
    if (trajin.SetupTrajRead(crdname, trajArgs, &parm, false)) return 1;
    TrajFrame.SetupFrameV(parm.Atoms(), trajin.TrajCoordInfo());
    trajin.BeginTraj(false);
    if (trajin.ReadTrajFrame(0, TrajFrame)) return 1;
    trajin.EndTraj();
  } else {
    // Assume Amber restart from STDIN
    // Check that input is from a redirect.
    if ( isatty(fileno(stdin)) ) {
      mprinterr("Error: No coordinates specified with '-c' and no STDIN '<'.\n");
      return 1;
    }
    Traj_AmberRestart restartIn;
    restartIn.SetDebug( debug );
    //restartIn.processReadArgs( trajArgs );
    int total_frames = restartIn.setupTrajin("", &parm);
    if (total_frames < 1) return 1;
    TrajFrame.SetupFrameV(parm.Atoms(), restartIn.CoordInfo());
    if (restartIn.openTrajin()) return 1;
    if (restartIn.readFrame(0, TrajFrame)) return 1;
    restartIn.closeTraj();
  }
  if (ctr_origin) 
    TrajFrame.CenterOnOrigin(false);
  // Output coords
  Trajout_Single trajout;
  trajArgs.SetList( aatm + bres + pqr + title, " " );
  if ( trajout.PrepareStdoutTrajWrite(trajArgs, &parm, fmt) ) return 1;
  trajout.WriteSingle(0, TrajFrame);
  trajout.EndTraj();
  return 0;
}
