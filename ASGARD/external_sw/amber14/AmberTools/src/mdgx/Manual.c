#include <stdio.h>
#include <string.h>
#include "Manual.h"

/***=======================================================================***/
/*** HorizontalRule: prints a horizontal bar to whatever file (including   ***/
/***                 stdout) specified.                                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   outp:   the target file                                             ***/
/***   int n:  the number of carriage returns to add after the rule        ***/
/***=======================================================================***/
void HorizontalRule(FILE *outp, int n)
{
  int i;

  fprintf(outp, "<++>---------------------------------------------------------"
	 "--------------<++>\n");
  for (i = 0; i < n; i++) {
    fprintf(outp, "\n");
  }
}

/***=======================================================================***/
/*** PrintSplash: print the splash lines for mdgx.  This is where to put   ***/
/***              information about the program's primary authorship and   ***/
/***              copyrights.                                              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   outp:  the target file                                              ***/
/***=======================================================================***/
void PrintSplash(FILE *outp)
{
  HorizontalRule(outp, 0);
  fprintf(outp,
          "<++> mdgx: A molecular dynamics engine in the AMBER suite of "
          "programs      <++>\n"
          "<++>                                                             "
          "          <++>\n"
          "<++> Written by David S. Cerutti, Case Group (2009)              "
          "          <++>\n");
  HorizontalRule(outp, 1);
}

/***=======================================================================***/
/*** PrintVADesc: this function prints a variable name, alias, and         ***/
/***              description using the specified formatting.              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   leadspace: the number of leading white space characters             ***/
/***   vname:     the variable name                                        ***/
/***   vnlen:     the amount of space to give the variable name            ***/
/***   valias:    the variable alias                                       ***/
/***   valen:     the amount of space to give the variable alias           ***/
/***   vdesc:     the variable description                                 ***/
/***   vdlen:     the amount of space to give the variable description (if ***/
/***                the description is longer than the amount of space     ***/
/***                alotted, additional lines will be used)                ***/
/***   vdindent:  the amount of indentation to apply to extra lines of the ***/
/***                variable description                                   ***/
/***   outp:      the output file                                          ***/
/***=======================================================================***/
void PrintVADesc(int leadspace, char* vname, int vnlen, char* valias,
		 int valen, char* vdesc, int vdlen, int vdindent, FILE *outp)
{
  int i, j, endfound, vpivot;
  char scratch[4096];

  /*** Print leading spaces ***/
  for (i = 0; i < leadspace; i++) {
    fprintf(outp, " ");
  }

  /*** Extend the variable name, then print it ***/
  j = strlen(vname);
  strcpy(scratch, vname);
  for (i = j; i < vnlen; i++) {
    scratch[i] = ' ';
  }
  scratch[vnlen] = '\0';
  fprintf(outp, "%s", scratch);

  /*** Extend the alias name, then print it ***/
  j = strlen(valias);
  strcpy(scratch, valias);
  for (i = j; i < valen; i++) {
    scratch[i] = ' ';
  }
  scratch[valen] = '\0';
  fprintf(outp, "%s", scratch);

  /*** Parse the description ***/
  endfound = 0;
  j = 0;
  while (endfound == 0) {
    for (i = 0; i < vdlen; i++) {
      if (vdesc[j+i] == ' ') {
	i++;
	vpivot = j+i;
      }
      if (vdesc[j+i] == '\0') {
	vpivot = j+i;
	endfound = 1;
	break;
      }
    }
    strncpy(scratch, &vdesc[j], vpivot-j);
    scratch[vpivot-j] = '\0';
    fprintf(outp, "%s\n", scratch);
    if (j == 0 && endfound == 0) {
      vdlen -= vdindent;
    }
    j = vpivot;
    if (endfound == 0) {
      for (i = 0; i < leadspace + vnlen + valen + vdindent; i++) {
	fprintf(outp, " ");
      }
    }
  }
}

/***=======================================================================***/
/*** PrintParagraph: prints a paragraph within a specified width.  No      ***/
/***                 leading white space or other columns are provided.    ***/
/***                 However, each paragraph is terminated by printing one ***/
/***                 additional carriage return (so that paragaphs are     ***/
/***                 separated by a full blank line).                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   vpar:   the paragraph string (carriage returns should not be        ***/
/***           included in the string)                                     ***/
/***   width:  the width of the text to print (the text will be broken     ***/
/***           up at whitespace characters to prevent lines from exceeding ***/
/***           this width)                                                 ***/
/***   outp:   the target file for this paragraph                          ***/
/***=======================================================================***/
void PrintParagraph(char *vpar, int width, FILE *outp)
{
  int i, j, vpivot, endfound;
  char scratch[4096];

  endfound = 0;
  j = 0;
  while (endfound == 0) {
    for (i = 0; i < width; i++) {
      if (vpar[j+i] == ' ') {
	i++;
        vpivot = j+i;
      }
      if (vpar[j+i] == '\0') {
        vpivot = j+i;
        endfound = 1;
        break;
      }
    }
    strncpy(scratch, &vpar[j], vpivot-j);
    scratch[vpivot-j] = '\0';
    fprintf(outp, "%s\n", scratch);
    j = vpivot;
  }
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** PrintUsage: print a brief set of usage instructions, to help users    ***/
/***             get started.                                              ***/
/***=======================================================================***/
void PrintUsage()
{
  printf("Usage:\n"
         ">> mdgx -i <input file>         (simplest operation)\n"
         ">> mdgx <arguments>             (operation similar to SANDER)\n"
         ">> mdgx <documentation name>    (prints documentation)\n\n");
  PrintParagraph("Command-line information may be entered using the same "
		 "arguments as the SANDER and PMEMD programs.  Alternatively, "
		 "all input may be provided in a single file using the -i "
		 "argument.  Any of the following arguments may also be used "
		 "to obtain additional documentation:", 79, stdout);
  printf("  -INPUT:   print all command line input options\n"
         "  -IFILE:   documentation on input file format\n"
	 "  -FILES:   print descriptions of &files namelist variables (these "
	 "may also\n"
	 "            be entered as command line input)\n"
         "  -CNTRL:   print descriptions of &cntrl namelist variables (most "
         "are similar\n"
         "            to SANDER variables, but some are unique to mdgx and "
         "some SANDER\n"
         "            variables are not supported)\n"
         "  -EWALD:   print &ewald namelist variables\n"
         "  -FORCE:   print &force namelist variables\n"
         "  -FITQ:    print &fitq (charge fitting) namelist variables\n"
         "  -PARAM:   print &param (bonded term fitting) namelist variables\n"
	 "  -IPOLQ:   print &ipolq (Implicitly Polarized Charge) namelist "
	 "variables\n"
         "  -ATTR:    attributions of certain aspects of the code, with "
	 "references\n\n");
}

/***=======================================================================***/
/*** PrintCommandLineInputOptions: this function essentially reproduces    ***/
/***                               what the AMBER manual already has for   ***/
/***                               SANDER and PMEMD command line input,    ***/
/***                               but since mdgx has some new features    ***/
/***                               it is necessary to have independent     ***/
/***                               documentation.                          ***/
/***=======================================================================***/
void PrintCommandLineInputOptions()
{
  PrintSplash(stdout);
  PrintVADesc(2, "-O", 5, " ", 2, "Overwrite output files if they exist "
	      "(appending files with the -A option found in SANDER is "
	      "currently not supported)", 71, 0, stdout);
  PrintVADesc(2, "-i", 5, " ", 2, "(input) control data for an energy "
	      "minimization / molecular dynamics run", 71, 0, stdout);
  PrintVADesc(2, "-o", 5, " ", 2, "(output) human-readable state information "
	      "and diagnostics", 71, 0, stdout);
  PrintVADesc(2, "-p", 5, " ", 2, "(input) molecular topology file (AMBER "
	      "prmtop format)", 71, 0, stdout);
  PrintVADesc(2, "-p2", 5, " ", 2, "(input) molecular topology file; if "
	      "thermodynamic integration is active, this topology describes "
	      "the final state while the original topology describes the "
	      "initial state", 71, 0, stdout);
  PrintVADesc(2, "-xpt", 5, " ", 2, "(input) extra points rule file directing "
	      "mdgx to add extra points to the topology at run time", 71, 0,
	      stdout);
  PrintVADesc(2, "-xpt2", 5, " ", 2, "(input) extra points rule file for the "
	      "topology specified by the -p2 flag", 71, 0, stdout);
  PrintVADesc(2, "-c", 5, " ", 2, "(input) initial coordinates (and, "
	      "optionally, velocities) and periodic box dimensions", 71, 0,
	      stdout);
  PrintVADesc(2, "-c2", 5, " ", 2, "(input) initial coordinates; if "
	      "thermodynamic integration is active, this second set of input "
	      "coordinates pertains to the initial coordinates of the system "
	      "in its final state as the mixing parameter lambda goes to 1",
	      71, 0, stdout);
  PrintVADesc(2, "-d", 5, " ", 2, "(output) comprehensive force / energy "
	      "report file", 71, 0, stdout);
  PrintVADesc(2, "-x", 5, " ", 2, "(output) coordinate trajectory file", 71,
	      0, stdout);
  PrintVADesc(2, "-x2", 5, " ", 2, "(output) coordinate trajectory file; only "
	      "used when thermodynamic integration is active", 71, 0, stdout);
  PrintVADesc(2, "-v", 5, " ", 2, "(output) velocity trajectory file", 71, 0,
	      stdout);
  PrintVADesc(2, "-v2", 5, " ", 2, "(output) velocity trajectory file; only "
	      "used when thermodynamic integration is active", 71, 0, stdout);
  PrintVADesc(2, "-e", 5, " ", 2, "(output) energy data over trajectory", 71,
	      0, stdout);
  PrintVADesc(2, "-r", 5, " ", 2, "(output) checkpoint (and final) "
	      "coordinates and periodic unit cell dimensions from energy "
	      "minimization runs, plus final velocities from molecular "
	      "dynamics runs\n", 71, 0, stdout);
  PrintVADesc(2, "-r2", 5, " ", 2, "(output) checkpoint file; only used when "
	      "thermodynamic integration is active\n", 71, 0, stdout);
}

/***=======================================================================***/
/*** PrintInputFormat: helpful documentation on the format of mdgx input   ***/
/***                   files.                                              ***/
/***=======================================================================***/
void PrintInputFormat()
{
  PrintSplash(stdout);
  PrintParagraph("The typical mdgx input file is designed to look very much "
		 "like a typical SANDER input file.  However, there are some "
		 "key differences implemented to make the mdgx input file "
		 "format more flexible and the control data more intuitive.",
		 79, stdout);
  PrintParagraph("The only significant restriction introduced to the mdgx "
		 "input file format is that different segments of the input "
		 "file must begin with the namelist identifier (e.g. &cntrl, "
		 "&ewald) on its own line, and be terminated with the "
		 "identifier &end, again on its own line.", 79, stdout);
  printf("Here is an example input file:\n\n"
	 "&files\n"
	 "  -p       Tests/ions.top\n"
	 "  -c       Tests/ions.min\n"
	 "  -rst     Tests/ionsMDGX\n"
	 "  -rstsuff .rst\n"
	 "&end\n\n"
	 "&cntrl\n"
	 "  DoRATTLE = 1,   LJ14Fac = 2.0,   Elec14Fac = 1.2,\n"
	 "  ElecCut = 9.0,  vdw_cutoff = 15.0,\n"
	 "  dt = 0.001,   nstlim 500000,  nfistep = 1000,\n"
	 "  ntpr = 1000,   ntwr 1000,  ntwf = 1000,\n"
	 "  Temperature = 0.0,\n"
	 "  SplnSpc = 0.015625,\n"
	 "&end\n\n"
	 "&ewald\n"
	 "  ordr1 = 4,\n"
	 "  ordr2 = 4,\n"
	 "  ordr3 = 4,\n"
	 "  nfft1 = 64,\n"
	 "  nfft2 = 64,\n"
	 "  nfft3 = 64;\n"
	 "&end\n\n");
  PrintParagraph("Note the presence of the familiar <namelist identifier> "
		 "<arguments> <&end> format, carried over from SANDER. "
		 "However, mdgx includes new namelists such as &files (which "
		 "allows the bulk of the command line information to be given "
		 "in the input file) and new variables such as file suffixes "
		 "(for specifying multiple files in a single run).  There are "
		 "also aliases for the familiar (though sometimes "
		 "unintelligible) SANDER namelist variables, to provide a "
		 "format that will accept SANDER input while also permitting "
		 "users to write less cryptic command files.", 79, stdout);
  PrintParagraph("Another important change to the mdgx file format is that = "
		 "signs are no longer required between a variable name and "
		 "the desired value.  In fact, commas are no longer required "
		 "either, though they are useful for separating different "
		 "attributes in each namelist.", 79, stdout);
}

/***=======================================================================***/
/*** PrintFilesNamelistVariables: this function provides documentation on  ***/
/***                              &files namelist variables, with whatever ***/
/***                              aliases may be available.                ***/
/***=======================================================================***/
void PrintFilesNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The files namelist is a means of specifying many input and "
		 "output files, which would otherwise be given on the command "
		 "line, in the input file along with other namelists. "
		 "However, if files are specified on the command line, they "
		 "will take precedence over any data in the &files namelist.",
		 79, stdout);
  PrintParagraph("The &files namelist can also be used to specify suffixes "
		 "for output files, in the event that multiple sequential "
		 "output files are to be generated in a single run.  This can "
		 "be useful for runs on managed resources that do not allow "
		 "multiple initializations of a program in a single job "
		 "submission, or for managing many segments of a very long "
		 "trajectory.  The suffixes only come into use if the "
		 "variable nfistep (alias FileStepCount) is set in the &cntrl "
		 "namelist.  Otherwise, only one file of each output type "
		 "will be written and only the base file names will be used. "
		 "File overwriting has slightly different behavior when "
		 "multiple sequential output files are specified.  If "
		 "overwriting is activated in such a case, mdgx will search "
		 "for the lowest file number such that a complete state "
		 "information file (specified by the -o variable) or a "
		 "complete restart file (specified by the -r variable) are "
		 "unavailable.  The dynamics will begin, or resume, at that "
		 "point.", 79, stdout);
  printf("  Name      Alias      Description\n"
	 " ------ ------------- ---------------------------------------------"
	 "-------------\n");
  PrintVADesc(1, "-p", 7, "Topology", 14, "(input) molecular topology file "
	      "(AMBER prmtop format)", 58, 2, stdout);
  PrintVADesc(1, "-p2", 7, "Topology2", 14, "(input) molecular topology file; "
	      "if thermodynamic integration is active, this topology "
	      "describes the final state while the original topology "
	      "describes the initial state", 58, 2, stdout);
  PrintVADesc(1, "-xpt", 7, "EPRules", 14, "(input) extra points rule file "
	      "directing mdgx to add extra points to the topology at run time",
	      58, 2, stdout);
  PrintVADesc(1, "-xpt2", 7, "EPRules2", 14, "(input) extra points rule file "
	      "for the topology specified by the -p2 flag", 58, 2, stdout);
  PrintVADesc(1, "-c", 7, "StartCrd", 14, "(input) initial coordinates (and, "
	      "optionally, velocities) and periodic unit cell size (mdgx does "
	      "not run with non-periodic unit cells)", 58, 2, stdout);
  PrintVADesc(1, "-d", 7, "ForceDump", 14, "(output) a comprehensive force "
	      "and energy report; this file is analogous to forcedump.dat as "
	      "produced by SANDER", 58, 2, stdout);
  PrintVADesc(1, "-rrp", 7, "ResReport", 14, "(output) a complete description "
	      "of the various residue types in the system; does not include "
	      "information on connections between residues", 58, 2, stdout);
  PrintVADesc(1, "-o", 7, "OutputBase", 14, "(output) human-readable state "
	      "information and diagnostics", 58, 2, stdout);
  PrintVADesc(1, "-e", 7, "EnergyBase", 14, "(output) energy data over "
	      "trajectory", 58, 2, stdout);
  PrintVADesc(1, "-x", 7, "CrdTrajBase", 14, "(output) coordinate trajectory "
	      "file; coordinate sets saved at specified time intervals", 58,
	      2, stdout);
  PrintVADesc(1, "-v", 7, "VelTrajBase", 14, "(output) velocity trajectory "
	      "file; velocity sets saved at specified time intervals", 58, 2,
	      stdout);
  PrintVADesc(1, "-f", 7, "FrcTrajBase", 14, "(output) force trajectory "
	      "file; forces on all particles saved at specified time"
	      "intervals", 58, 2, stdout);
  PrintVADesc(1, "-r", 7, "RestartBase", 14, "(output) checkpoint (and final) "
	      "coordinates and periodic unit cell dimensions from energy "
	      "minimization runs, plus final velocities from molecular"
	      "dynamics runs", 58, 2, stdout);
  PrintVADesc(1, "-osf", 7, "OutputSuff", 14, "output state information data "
	      "file suffix (if this or other suffixes are not specified, the "
	      "base name is taken to be the complete file name)", 58, 2,
	      stdout);
  PrintVADesc(1, "-esf", 7, "EnergySuff", 14, "Energy data file suffix", 58,
	      2, stdout);
  PrintVADesc(1, "-xsf", 7, "CrdTrajSuff", 14, "Coordinate trajectory file "
	      "suffix", 58, 2, stdout);
  PrintVADesc(1, "-vsf", 7, "VelTrajSuff", 14, "Velocity trajectory file "
	      "suffix", 58, 2, stdout);
  PrintVADesc(1, "-fsf", 7, "FrcTrajSuff", 14, "Force trajectory file suffix",
	      58, 2, stdout);
  PrintVADesc(1, "-rsf", 7, "RestartSuff", 14, "Restart file suffix\n", 58, 2,
	      stdout);
}

/***=======================================================================***/
/*** PrintCntrlNamelistVariables: this function provides documentation on  ***/
/***                              &cntrl namelist variables, with whatever ***/
/***                              aliases may be available.                ***/
/***=======================================================================***/
void PrintCntrlNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The &cntrl namelist is required for any SANDER run, and is "
		 "present with little modification in mdgx.  Most SANDER "
		 "input files can therefore be read directly by mdgx or "
		 "adapted without much effort.  However, there are some "
		 "variables that have either been replaced or discarded.  The "
		 "following list describes all variables available in the "
		 "mdgx &cntrl namelist.", 79, stdout);
  PrintParagraph("Molecular Dynamics Variables:", 79, stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- -----------------------------------------"
         "-------------\n");
  PrintVADesc(1, "imin", 11, "RunMode", 14, "The run mode (0 = molecular "
	      "dynamics, 1 = energy minimization, 2 = force computation for "
	      "the input coordinates)", 54, 2, stdout);
  PrintVADesc(1, "irest", 11, "RestartMD", 14, "Set to 1 to request that "
	      "molecular dynamics be restarted using velocities given in the "
	      "input coordiantes file; set to 0 to request that initial " 
	      "velocities be assigned to random values in a Maxwell "
	      "distribution", 54, 2, stdout);
  PrintVADesc(1, "ioutfm", 11, "CoordFormat", 14, "The format of trajectory "
              "output coordinates, 0 (default) being ascii format with three "
	      "decimal places and 1 being binary NetCDF format, specifying "
	      "all coordinates to single floating point precision", 54, 2,
	      stdout);
  PrintVADesc(1, "nstlim", 11, "StepCount", 14, "Number of MD steps to be "
	      "performed", 54, 2, stdout);
  PrintVADesc(1, "nfistep", 11, "FileStepCount", 14, "Length of each segment "
              "of the trajectory; the number of segments is nstlim / nfistep",
	      54, 2, stdout);
  PrintVADesc(1, "nscm", 11, "ZeroMomentum", 14, "Removal of translational "
	      "and rotational center-of-mass (COM) motion at regular "
	      "intervals", 54, 2, stdout);
  PrintVADesc(1, "t", 11, "StartTime", 14, "Time at the start of the "
	      "simulation (picoseconds); this parameter is for user reference "
	      "and otherwise only affects the time entered in the diagnostics "
	      "files", 54, 2, stdout);
  PrintVADesc(1, "dt", 11, "TimeStep", 14, "Time step (picoseconds); the "
	      "recommended MAXIMUM is 0.0025 if bonds to hydrogen are "
	      "constrained, or 0.001 otherwise", 54, 2, stdout);
  PrintVADesc(1, "temp0", 11, "Temperature", 14, "Target temperature for a "
              "constant temperature simulation.  The default is 298.0.  This "
	      "value can be used to initialize velocities if a specific "
	      "initial temperature is not set.", 54, 2, stdout);
  PrintVADesc(1, "tempi", 11, "StartTemp", 14, "Initial temperature for a "
              "simulation.  The default is -100.0, which commands mdgx to use "
	      "temp0 as the initial temperature for things such as velocity "
	      "initialization.  If a positive value of tempi is specified, "
	      "tempi will be used to initialize velocities.", 54, 2, stdout);
  PrintVADesc(1, "ntt", 11, "Thermostat", 14, "Thermostat specification.  "
	      "Numerical values of 0 (default, no thermostat), 1 (Berendsen "
	      "thermostat), 2 (Andersen thermostat), 3 (Langevin integrator) "
	      "and 4 (Nose-Hoover thermostat) are supported.", 54, 2, stdout);
  PrintVADesc(1, "ntp", 11, "CoordRscl", 14, "Coordinate rescaling "
	      "specification.  Numerical values of 0 (default, no rescaling, "
	      "constant volume), 1 (isotropic rescaling), and 2 "
	      "(anisotropic rescaling) are supported.", 54, 2, stdout);
  PrintVADesc(1, "barostat", 11, "Barostat", 14, "Barostat style.  Numerical "
	      "values of 1 (Berendsen, default), and 2 (Monte-Carlo) are "
	      "supported.", 54, 2, stdout);
  PrintVADesc(1, "pres0", 11, "Pressure", 14, "Target pressure for a constant "
              "pressure simulation.  The default is 1.0 bar.", 54, 2, stdout);
  PrintVADesc(1, "tautp", 11, "BerendsenTC", 14, "Time constant for "
	      "Berendsen temperature coupling.  Default value is 0.4 ps.", 54,
              2, stdout);
  PrintVADesc(1, "gamma_ln", 11, "LangevinFreq", 14, "Langevin collision "
	      "frequency in events / ps, when ntt=3.  Default is 0.", 54, 2,
	      stdout);
  PrintVADesc(1, "taup", 11, "BerendsenPC", 14, "Compressibility constant for "
	      "Berendsen pressure coupling.  Default value is 4.4e-5 / bar.",
	      54, 2, stdout);
  PrintVADesc(1, "tauthv", 11, "HooverTC", 14, "Time constant for Hoover "
	      "temperature coupling.  Default value is 1.0 ps.", 54, 2,
	      stdout);
  PrintVADesc(1, "tauphv", 11, "HooverPC", 14, "Compressibility constant for "
              "Hoover pressure coupling.  Default value is 1.0 / bar.", 54, 2,
	      stdout);
  PrintVADesc(1, "mccomp", 11, "MCBarostatPC", 14, "Coordinate rescaling "
	      "factor for isotropic Monte-Carlo pressure coupling.  Default "
	      "is 2.0e-3, to rescale the volume by +/- 1/10th of 1%.", 54, 2,
	      stdout);
  PrintVADesc(1, "mccompx", 11, "MCBarostatPCX", 14, "Coordinate rescaling "
	      "factor for anisotropic Monte-Carlo pressure coupling in the X "
	      "direction.  Default is 2.0e-3, to rescale the volume by +/- "
	      "1/10th of 1%.", 54, 2, stdout);
  PrintVADesc(1, "mccompy", 11, "MCBarostatPCY", 14, "Coordinate rescaling "
	      "factor for anisotropic Monte-Carlo pressure coupling in the Y "
	      "direction.  Default is to perform only isotropic rescaling "
	      "based on the factor for the X direction.", 54, 2, stdout);
  PrintVADesc(1, "mccompz", 11, "MCBarostatPCZ", 14, "Coordinate rescaling "
	      "factor for anisotropic Monte-Carlo pressure coupling in the Z "
	      "direction.  Default is to perform only isotropic rescaling "
	      "based on the factor for the X direction.", 54, 2, stdout);
  PrintVADesc(1, "mcbfrq", 11, "MCBarostatFrq", 14, "Step frequency for "
	      "applying the Monte-Carlo barostat, if this barostat is "
	      "activated.  Default 100.", 54, 2, stdout);
  PrintVADesc(1, "vrand", 11, "RandomReset", 14, "Time constant for Andersen "
	      "temperature coupling.  This is specified as an integer "
	      "denoting the number of time steps, and as such is related to "
	      "the time step size dt.  Default value is 1000.", 54, 2, stdout);
  PrintVADesc(1, "ig", 11, "RandomSeed", 14, "The random seed for velocity "
	      "initialization and thermostats which may require it", 54, 2,
	      stdout);
  PrintVADesc(1, "es_cutoff", 11, "ElecCut", 14, "The electrostatic direct "
	      "space cutoff (Angstroms)", 54, 2, stdout);
  PrintVADesc(1, "vdw_cutoff", 11, "VdwCut", 14, "The van-der Waals (direct "
	      "space) cutoff (Angstroms)", 54, 2, stdout);
  PrintVADesc(1, "cut", 11, "DirectCutoff", 14, "The general (van-der Waals "
              "and electrostatic) direct space cutoff (Angstroms).  This "
	      "value, if indicated, will override es_cutoff and vdw_cutoff.",
	      54, 2, stdout);
  PrintVADesc(1, "rigidbond", 11, "DoRATTLE", 14, "Set to 1 to activate "
	      "RATTLE bond length constraints", 54, 2, stdout);
  PrintVADesc(1, "rigidwat", 11, "DoSETTLE", 14, "Set to 1 to activate "
	      "SETTLE water geometry constraints", 54, 2, stdout);
  PrintVADesc(1, "tol", 11, "RattleTol", 14, "Tolerance for RATTLE bond "
	      "length constraints", 54, 2, stdout);
  PrintVADesc(1, "scee", 11, "Elec14Fac", 14, "The electrostatic 1-4 "
	      "interaction scaling factor", 54, 2, stdout);
  PrintVADesc(1, "scnb", 11, "Vdw14Fac", 14, "The van-der Waals 1-4 "
	      "interaction scaling factor", 54, 2, stdout);
  printf("\n");
  PrintParagraph("Thermodynamic integration control variables:", 79, stdout);
  printf("  Name        Alias        Description\n"
         " ------ ----------------- -----------------------------------------"
         "-------------\n");
  PrintVADesc(1, "icfe", 11, "RunTI", 14, "Flag to activate thermodynamic "
	      "integration.  Default 0 (no TI); set to 1 for active.", 54, 2,
	      stdout);
  PrintVADesc(1, "clambda", 11, "MixFactor", 14, "The mixing parameter, "
	      "(1-L)^k of state 1 and 1-(1-L)^k of state 2", 54, 2, stdout);
  PrintVADesc(1, "klambda", 11, "MixOrder", 14, "The order of the mixing "
	      "parameter, (1-L)^k of state 1 and 1-(1-L)^k of state 2", 54, 2,
	      stdout);
  PrintVADesc(1, "nsynch", 11, "SynchTI", 14, "Explicit synchronization of "
	      "trajectory coordinates will occur every nsynch steps (default "
	      "1000).  Set nsynch to 0 to disable this feature.", 54, 2,
	      stdout);
  printf("\n");
  PrintParagraph("Output control variables:", 79, stdout);
  printf("  Name        Alias        Description\n"
         " ------ ----------------- -----------------------------------------"
         "-------------\n");
  PrintVADesc(1, "ntpr", 7, "WriteDiagnostics", 18, "Diagnostics and state "
	      "information output frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwx", 7, "WriteCoordinates", 18, "Trajectory coordinates "
	      "will be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwv", 7, "WriteCoordinates", 18, "Trajectory velocities "
	      "will be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwf", 7, "WriteCoordinates", 18, "Trajectory forces will "
	      "be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwr", 7, "WriteCoordinates", 18, "Trajectory forces will "
	      "be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "tchk", 7, "TopologyCheck", 18, "Active by default, set to 0 "
	      "to turn off topology and conformation checking at the start of "
	      "each run segment", 54, 2, stdout);
}

/***=======================================================================***/
/*** PrintEwaldNamelistVariables: this function provides documentation on  ***/
/***                              &ewald namelist variables, with whatever ***/
/***                              aliases may be available.                ***/
/***=======================================================================***/
void PrintEwaldNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("All variables in the &ewald namelist have default values, "
		 "so use of this namelist is optional, but optimization of "
		 "the parameters in this section can be very helpful for "
		 "performing the most efficient molecular simulations.  The "
		 "SANDER and PMEMD programs warn users not to modify "
		 "variables in this section without significant experience; "
		 "what is most important is a clear understanding of all the "
		 "variables and how they will affect the accuracy of "
		 "computed forces.  The state information / diagnostics "
		 "output file will print verbose explanations of the "
		 "consequences of any variables that are changed, so "
		 "modification of these variables, with careful reading of "
		 "the output and checks on the accuracy of computed forces, "
		 "should be safe.", 79, stdout);
  PrintParagraph("Smooth Particle Mesh Ewald (SPME) Variables:", 79, stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- ------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "dsum_tol", 11, "DSumTol", 14, "The direct sum tolerance; "
	      "at the direct space cutoff es_cutoff (see the &cntrl namelist)"
	      ", the ratio between the interaction energy of two point "
	      "charges and two Gaussian smeared charges of the same magnitude "
	      "will differ from 1 by dsum_tol; this variable controls the "
	      "accuracy of the electrostatic direct space sum", 54, 2, stdout);
  PrintVADesc(1, "sigma", 11, "Sigma", 14, "The width (root mean squared "
	      "deviation) of spherical Gaussians used to smooth out the "
	      "charge density; sigma = ew_coeff / 2; if this value is "
	      "supplied by the user it will supercede any entry for ew_coeff "
	      "and dsum_tol will be calculated from sigma; otherwise, sigma "
	      "and ew_coeff will be calculated from dsum_tol", 54, 2, stdout);
  PrintVADesc(1, "ew_coeff", 11, "EwaldCof", 14, "The Ewald coefficient; this "
	      "quantity has a name only because it appears frequently in "
	      "the Smooth Particle Mesh Ewald formulas; physically it makes "
	      "more sense to consider the Gaussian charge width sigma", 54, 2,
	      stdout);
  PrintVADesc(1, "eetbdns", 11, "SplnSpc", 14, "The discretization of the "
	      "erfc(beta*r)/r force and energy spline computation tables", 54,
	      2, stdout);
  PrintVADesc(1, "rho", 11, "MaxDensity", 14, "The maximum expected density "
	      "of the system, g/mL.  Default 2.0, increase to raise the "
	      "maximum storage in the direct-space decomposition cell grid.",
	      54, 2, stdout);
  PrintVADesc(1, "nfft[123]", 11, "MeshDim[XYZ]", 14, "The number of mesh "
	      "points in the X, Y, or Z dimensions, respectively (if the unit "
	      "cell / simulaton box is orthorhombic), or otherwise the number "
	      "of mesh points along the 1st, 2nd, and 3rd unit cell vectors",
	      54, 2, stdout);
  PrintVADesc(1, "ordr[123]", 11, "Order[XYZ]", 14, "The order of particle "
	      "<--> mesh interpolation along the 1st, 2nd, and 3rd unit cell "
	      "vectors", 54, 2, stdout);
  PrintVADesc(1, "order", 11, "Order", 14, "Sets the interpolation order "
	      "isotropically along all unit cell vectors to the specified "
	      "value (supercedes ordr[123])", 54, 2, stdout);
  printf("\n");
  PrintParagraph("Long-ranged van-der Waals control parameters:", 79, stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- ------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "vdwmeth", 11, "vdwMethod", 14, "The method for computing "
	      "long-ranged van der Waals interactions.  The default of 1 "
	      "implies the inclusion of a homogeneity assumption in the "
	      "long-ranged component of the van-der Waals interactions; this "
	      "correction changes the computed energy and pressure, but not "
	      "forces, and would therefore not affect dynamics in a "
	      "simulation at constant volume. A value of zero removes any "
	      "such correction and make the van-der Waals energy depend "
	      "solely on the pairwise interactions.", 54, 2, stdout);
  printf("\n");
  PrintParagraph("Multi-Level Ewald (MLE) Variables:", 79, stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- ------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "nlev", 11, "EwaldLevels", 14, "The number of levels in the "
	      "multi-level mesh hierarchy (standard Smooth Particle Mesh "
	      "Ewald has one level); setting this variable to any value "
	      "greater than one will activate MLE (maximum 4)", 54, 2, stdout);
  PrintVADesc(1, "lpad[123]", 11, "Padding[123]", 14, "The number of layers "
	      "of \"padding\" for meshes at levels 1, 2, and 3; the "
	      "reciprocal space pair potential will be represented exactly on "
	      "mesh level 1 up to (and including) lpad1 grid points from the "
	      "source, and will be represented to different degrees of "
	      "resolution (see cfac[234]) on higher mesh levels, up to "
	      "(and including) lpad1 + lpad2 or lpad1 + lpad2 + lpad3 layers "
	      "from the source (note that the highest mesh level is not "
	      "padded as it involves only one convolution over the entire "
	      "simulation box); higher values of lpad[123] will produce more "
	      "accurate results (see also ggordr)", 54, 2, stdout);
  PrintVADesc(1, "cfac[234]", 11, "Spread[234]", 14, "The coarsening factor "
	      "for higher mesh levels; by definition, cfac1 is 1; generally, "
	      "it is advisable to set cfac2 to 1.5-2.0, and to set cfac3 or "
	      "cfac4 (if even higher mesh levels are in use) to increasingly "
	      "large numbers; although cfac is a real number, but it must be "
	      "specified such that the mesh size is an integer multiple of "
	      "cfac in every dimension", 54, 2, stdout);
  PrintVADesc(1, "ggordr", 11, "GridOrder", 14, "The order of grid <--> grid "
	      "B-Spline interpolation; higher orders of interpolation will "
	      "produce more accurate results for given values of lpad[123]",
	      54, 2, stdout);
}

/***=======================================================================***/
/*** PrintForceNamelistVariables: the force namelist was added to support  ***/
/***                              customization of detailed force report   ***/
/***                              files.                                   ***/
/***=======================================================================***/
void PrintForceNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("By default, all available information is printed to the "
		 "force report file; all of the variables listed below are "
		 "set to 1 by default.  Most of this information is not "
		 "needed, however, so much of the output can be suppressed by "
		 "setting these variables to 0.", 79, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "var", 11, "VarName", 14, "Information is dumped into a "
	      "Matlab-readable force report file; this specifies the name of "
	      "the variable that Matlab will use to store the force "
	      "information if it reads the report; it may be useful to read "
	      "multiple reports at once, and compare the results using Matlab",
	      54, 2, stdout);
  PrintVADesc(1, "dumpcrd", 11, "DumpCoord", 14, "Flag to dump the "
	      "coordinates of the system", 54, 2, stdout);
  PrintVADesc(1, "dumpbond", 11, "DumpBond", 14, "Flag to dump the forces due "
	      "to bonded (1-2) interactions", 54, 2, stdout);
  PrintVADesc(1, "dumpangl", 11, "DumpAngl", 14, "Flag to dump the forces due "
	      "to angle interactions", 54, 2, stdout);
  PrintVADesc(1, "dumpdihe", 11, "DumpDihe", 14, "Flag to dump the dihedral "
	      "forces", 54, 2, stdout);
  PrintVADesc(1, "dumpdelec", 11, "DumpDElec", 14, "Flag to dump direct sum "
	      "electrostatic forces", 54, 2, stdout);
  PrintVADesc(1, "dumprelec", 11, "DumpRElec", 14, "Flag to dump reciprocal "
	      "sum electrostatic forces", 54, 2, stdout);
  PrintVADesc(1, "dumpvdw", 11, "DumpVdw", 14, "Flag to dump van-der Waals "
	      "forces", 54, 2, stdout);
  PrintVADesc(1, "dumpall", 11, "DumpAll", 14, "Flag to dump total (summed) "
	      "forces", 54, 2, stdout);
}

/***=======================================================================***/
/*** PrintFitqNamelistVariables: the fitq namelist was added to support    ***/
/***                             greatly expanded extra points features.   ***/
/***=======================================================================***/
void PrintFitqNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The basic concept of fitting charges to reproduce the "
		 "electrostatic potential of a molecule, by finding the "
		 "solution with least squared error, in the presence of "
		 "restraints, is carried over from the original Kollmann "
		 "RESP.  This namelist provides tools for fitting charges in "
		 "a comprehensible manner, with exceptional user control over "
		 "the range of the fitting data.  Because mdgx is unique "
		 "among the current AMBER molecular dynamics engines for its "
		 "ability to use certain types of extra points (virtual "
		 "sites), this namelist is also useful for fitting new charge "
		 "models to accelerate parameter development in mdgx.", 79,
		 stdout);
  PrintParagraph("Many of these variables can be specified more than once, "
		 "and all instances will accumulate in the result.", 79,
		 stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "resp", 11, "RespPhi", 14, "File names and numerical weight "
	      "of an electrostatic potential to use in fitting.  The format "
	      "is <string1> <string2> <real1>, where string1 is a Gaussian "
	      "cubegen format file specifying the electrostatic potential "
	      "and molecular coordinates and Z-numbers appropriate to the "
	      "topology specified by string2 and real1 is the numerical "
	      "weight to be given to this conformation in the fit.  This "
	      "variable may be specified for all molecular conformations "
	      "and systems to be used in the fit.\n", 54, 2, stdout);
  PrintVADesc(1, "ipolq", 11, "IPolQPhi", 14, "File names and numerical "
	      "weight of a pair of electrostatic potentials to use in IPolQ "
	      "fitting.  The format is <string1> <string2> <string3> <real1>, "
	      "where string1 and string2 are Gaussian cubegen format files "
	      "relating to the system in vacuum and in a condensed-phase "
	      "environment (see the &ipolq namelist).  Note that the "
	      "molecular coordinates in both cubegen files must match.  As in "
	      "the resp variable format, string3 is the topology, and real1 "
	      "is the numerical weight of this conformation.  This extended "
	      "REsP procedure supports development of fixed-charge force "
	      "fields if one posits that the correct charges of a "
	      "non-polarizable model would sit halfway between the charges of "
	      "a fully polarized molecule in some solvent reaction field and "
	      "the charges of an unpolarized molecule in the gas phase.  This "
              "variable may be specified for all molecular conformations "
              "and systems to be used in the fit.\n", 54, 2, stdout);
  PrintVADesc(1, "conf", 11, "ConfFile", 14, "If specified, mdgx, will output "
	      "the first molecular conformation, complete with any added "
	      "virtual sites, in PDB format for inspection.  This is useful "
	      "for understanding exactly what model is being fitted.", 54, 2,
	      stdout);
  PrintVADesc(1, "eprules", 11, "EPRules", 14, "If specified, mdgx, will "
	      "output all fitted charges in the form of a Virtual Sites rule "
	      "file, which can be given as input to subsequent simulations to "
	      "modify the original prmtp and apply the fitted charge model.",
	      54, 2, stdout);
  PrintVADesc(1, "hist", 11, "HistFile", 14, "If specified, mdgx will print "
	      "a histogram of the distance of all points from the nearest "
	      "atom of the molecule.", 54, 2, stdout);
  PrintVADesc(1, "minq", 11, "MinimizeQ", 14, "Restrain the charges of a "
	      "group of atoms to zero by the weight given in minqwt.  This "
	      "variable may be specified many times for different groups.  "
	      "The groups are specified in ambmask format.", 54, 2, stdout);
  PrintVADesc(1, "equalq", 11, "EqualizeQ", 14, "Restrain the charges of a "
	      "group of atoms to have the same values.  Groups are specified "
	      "in ambmask format.  This variable may be repeatedly specified.",
	      54, 2, stdout);
  PrintVADesc(1, "sumq", 11, "SumQ", 14, "Restrain the total charge of a "
	      "group of atoms to have a particular value.  The group is "
	      "specified in ambmask format, followed by a real number "
	      "indicating the target total charge.  A stiff harmonic "
	      "restraint penalty is applied.  This variable may be repeatedly "
	      "specified.", 54, 2, stdout);
  PrintVADesc(1, "snap", 11, "MaxSnap", 14, "Charges will be adjusted to "
	      "eliminate miniscule deviations in the total charge of "
	      "specified groups.  This variable sets the maximum adjustment "
	      "that can be made, in increments of 1.0e-5 proton charges.", 54,
	      2, stdout);
  PrintVADesc(1, "tether", 11, "Tether", 14, "Restrain fitted charges not "
	      "zero (small values), but to their values in the input "
	      "topologies, to use the original force field as a guide.", 54, 2,
	      stdout);
  PrintVADesc(1, "tetherwt", 11, "TetherWeight", 14, "Weight used for "
	      "harmonic restraint of fitted charges to values given in the "
	      "original force field (in standard REsP), or for restraining "
	      "IPolQ charges to their vacuum counterparts (in IPolQ fitting).",
	      54, 2, stdout);
  PrintVADesc(1, "minqwt", 11, "MinQWeight", 14, "Weight used for restraining "
	      "values of charges to zero; as more and more fitting data is "
	      "included (either through a higher sampling density of the "
	      "electrostatic potential due to each molecular conformation or "
	      "additional molecular conformations) higher values of minqwt "
	      "may be needed to keep the fitted charges small.  However, with "
	      "more data the need to restrain charges may diminish as well.",
	      54, 2, stdout);
  PrintVADesc(1, "nfpt", 11, "FitPoints", 14, "The number of fitting points "
	      "to select from each electrostatic potential grid.  The points "
	      "nearest the molecule, which satisfy the limits set by the "
	      "solvent probe and point-to-point distances as defined below, "
	      "will be selected for the fit.  Default 1000.", 54, 2, stdout);
  PrintVADesc(1, "psig", 11, "ProbeSig", 14, "The Lennard-Jones sigma "
	      "parameter of the solvent probe.  Default 3.16435 (TIP4P "
	      "oxygen).", 54, 2, stdout);
  PrintVADesc(1, "peps", 11, "ProbeEps", 14, "The Lennard-Jones parameter of "
	      "the solvent probe.  Default 0.16275 (TIP4P oxygen).", 54, 2,
	      stdout);
  PrintVADesc(1, "parm", 11, "ProbeArm", 14, "The probe arm; points on the "
	      "electrostatic potential grid that would be inaccessible to the "
	      "solvent probe may still be included in the fit if they are "
	      "within the probe arm's reach.", 54, 2, stdout);
  PrintVADesc(1, "pnrg", 11, "StericLimit", 14, "The maximum Lennard-Jones "
	      "energy of the solvent probe at which a point will qualify for "
	      "inclusion in the fit.  Default 3.0 kcal/mol.", 54, 2, stdout);
  PrintVADesc(1, "maxmem", 11, "MaxMemory", 14, "The amount of memory "
	      "available to this job.  The test is actually the amount "
	      "allocated to the fitting and testing matrices, which is the "
	      "number of fitting data points times the number of independent "
	      "charges being fitted times 16 bytes.  The actual memory usage "
	      "will be slightly higher, but the matrices comprise the bulk of "
	      "it.\n", 54, 2, stdout);
  PrintVADesc(1, "flim", 11, "Proximity", 14, "The minimum proximity of any "
	      "two points to be included in the fit.  Default 0.4A.", 54, 2,
	      stdout);
  PrintVADesc(1, "hbin", 11, "HistogramBin", 14, "If hist is specified, mdgx "
	      "will print a histogram reporting the number of fitting points "
	      "falling within any particular distance of some atom of the "
	      "molecule.  This parameter controls the discretization of the "
	      "histogram.", 54, 2, stdout);
  PrintVADesc(1, "verbose", 11, "Verbose", 14, "Print information relating to "
	      "progress on the fitting run.  Default is to print such data.  "
	      "Set to zero to suppress output.", 54, 2, stdout);
}

/***=======================================================================***/
/*** PrintParamNamelistVariables: the param namelist supports development  ***/
/***                              of bonded terms, particularly torsion    ***/
/***                              potentials, with consideration to each   ***/
/***                              parameter possibly having multiple roles ***/
/***                              in many different systems.               ***/
/***=======================================================================***/
void PrintParamNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("Bonded terms can be fitted by a linear least squares "
		 "approach, tuning the stiffnesses of various spring-like "
		 "interactions with quantum data as the guide.  The fitting "
		 "functions are fundamentally the same as those used by the "
		 "charge fitting module, but the &param namelist can "
		 "accommodate many more systems.  The basic requirements of "
		 "this module are an Amber parameter file (i.e. parm99.dat), "
		 "an optional frcmod file (i.e. frcmod.ff99SB), and a list of "
		 "system topologies and coordinates coupled to energies "
		 "obtained by a consistent quantum mechanics treatment.  mdgx "
		 "will correlate all adjustable bonded terms found in any of "
		 "the systems by referencing the parameter files, then "
		 "compute a fitting matrix with one adjustable term per "
		 "column.  The fitted parameters are printed to a formatted "
		 "Amber parameter file, while statistics from the run, "
		 "including accuracy in each system and a comprehensive "
		 "analysis of the fitting data, are printed to mdout.", 79,
		 stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "sys", 11, "System", 14, "A fitting data point. This keyword "
              "must be followed by three items: the name of a topology file, "
              "the name of a corresponding coordinate file, and the energy of "
              "this system in the stated conformation.", 54, 2, stdout);
  PrintVADesc(1, "bonds", 11, "FitBonds", 14, "Requests a linear "
	      "least-squares fit for bond stiffnesses in the system.", 54, 2,
	      stdout);
  PrintVADesc(1, "angles", 11, "FitAngles", 14, "Requests a linear "
	      "least-squares fit for angle stiffnesses in the system.", 54, 2,
	      stdout);
  PrintVADesc(1, "torsions", 11, "FitTorsions", 14, "Requests a linear "
	      "least-squares fit for torsion stiffnesses in the system.", 54,
	      2, stdout);
  PrintVADesc(1, "fith", 11, "FitH", 14, "Request that a specific torsion "
	      "parameter be included in linear least-squares fitting.", 54, 2,
	      stdout);
  PrintVADesc(1, "fitscnb", 11, "FitLJ14", 14, "Requests a linear "
	      "least-squares fit for Lennard-Jones 1:4 scaling factors.",
	      54, 2, stdout);
  PrintVADesc(1, "fitscee", 11, "FitEE14", 14, "Requests a linear "
	      "least-squares fit for electrostatic 1:4 scaling factors.",
	      54, 2, stdout);
  PrintVADesc(1, "repall", 11, "ReportAll", 14, "Flag to activate output "
	      "of all parameters encountered during the fitting procedure, "
	      "including those that were not adjusted by the fit but "
	      "nonetheless contributed to the molecular mechanics energies.  "
	      "Default is 1 (write all parameters to the Amber parameter "
	      "file), appropriate for creating a parm##.dat file to specify a "
	      "new force field.  Set to 0 to create files more akin to "
	      "frcmod files.  If branching atom types, setting this input to "
	      "2 will print all parameters, fitted or not, included in energy"
              "calculations or not, in the output files.", 54, 2, stdout);
  PrintVADesc(1, "zeromm", 11, "FittedMMOnly", 14, "Flag to ignore any "
	      "contributions from unfitted energy terms during the fit.  "
	      "Useful for making force field adjustments.", 54, 2, stdout);
  PrintVADesc(1, "verbose", 11, "ShowProgress", 14, "Alert the user as to the "
	      "progress of the fitting procedure.  Runs involving thousands "
	      "of molecular conformations and hundreds of parameters can "
	      "generally be completed in a few minutes.  Default is 1 (ON).  "
	      "Set to zero to suppress output.", 54, 2, stdout);
  PrintVADesc(1, "elimsig", 11, "ElimOutliers", 14, "Flag to activate removal "
	      "of molecular conformations whose energies are far outside the "
	      "norm for other conformations of the same system.  Default 0 "
	      "(do not remove outliers).", 54, 2, stdout);
  PrintVADesc(1, "esigtol", 11 ,"EOutlier", 14, "Tolerance for deviation from "
	      "the mean energy value, specified as a function of the standard "
	      "deviation for all conformations of the same system.  "
	      "Conformations of a system which exceed this threshold will be "
	      "reported if verbose is set to 1, and removed from "
	      "consideration if elimsig is set to 1.  Default 5.0 sigmas.\n",
	      54, 2, stdout);
  PrintVADesc(1, "ctol", 11, "ConfTol", 14, "Energetic tolerance for strained "
	      "bonds or angles in the system.  If any bonded energy term "
	      "exceeds this tolerance, mdgx will assume that the atoms have "
	      "been entered in the wrong order according to the topology, and "
	      "try to rearrange the atoms in order to relax the system.", 54,
	      2, stdout);  
  PrintVADesc(1, "eunits", 11, "EnergyUnits", 14, "Units of the target energy "
	      "values.  Default kcal/mol.  Acceptable values include Hartree/"
	      "Atomic, kJ/kilojoules, and j/joules.  Case insensitive.", 54,
	      2, stdout);
  PrintVADesc(1, "accrep", 11, "AccReport", 14, "Accuracy report on the fit.  "
	      "Contains extensive analysis on the resulting parameters, "
	      "in MatLab format.", 54, 2, stdout);
  PrintVADesc(1, "title", 11, "ParmTitle", 14, "Parameter file title.  This "
	      "is not a file name, but rather the title appearing on the "
	      "first line of the printed file.", 54, 2, stdout);
  PrintVADesc(1, "scnb", 11, "Vdw14Fac", 14, "Sets a universal 1:4 scaling "
	      "factor for van-der Waals interactions.  Use this input to "
	      "change the scaling on all systems simultaneously.", 54, 2,
	      stdout);
  PrintVADesc(1, "scee", 11, "Elec14Fac", 14, "Sets a universal 1:4 scaling "
	      "factor for electrostatic interactions.  Use this input to "
	      "change the scaling on all systems simultaneously.", 54, 2,
	      stdout);
  PrintVADesc(1, "brst", 11 ,"BondRest", 14, "General value for harmonic "
	      "restraints on bond stiffness constants.", 54, 2, stdout);
  PrintVADesc(1, "arst", 11 ,"AngleRest", 14, "General value for harmonic "
	      "restraints on angle stiffness constants.", 54, 2, stdout);
  PrintVADesc(1, "hrst", 11 ,"DihedralRest", 14, "General value for harmonic "
	      "restraints on torsion stiffness constants.", 54, 2, stdout);
}

/***=======================================================================***/
/*** PrintIPolQNamelistVariables: the ipolq namelist expedites the         ***/
/***                              otherwise laborious process of computing ***/
/***                              a solvent reaction field potential for a ***/
/***                              molecule and then setting up quantum     ***/
/***                              calculations in the vacuum and condensed ***/
/***                              phases.                                  ***/
/***=======================================================================***/
void PrintIPolQNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The Implicitly Polarized Charge method is laborious to "
		 "implement manually, but this module and associated namelist "
		 "automate much of the process, prevent errors, and make it "
		 "easy to determine convergence in the calculations.  A "
		 "typical REsP procedure requires a series of conformations "
		 "of the molecule of interest; the IPolQ procedure requires "
		 "a series of conformations of the molecule in a condensed "
		 "phase environment, such as water.  With the automation "
		 "afforded by this namelist, users can derive IPolQ charges "
		 "nearly as easily as standard REsP charges.\n", 79, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "solute", 11, "SoluteMol", 14, "The solute molecule, "
	      "specified by an ambmask string.  This is the molecule of "
	      "interest for charge fitting, and will be immobilized during "
	      "the simulation.  This must be specified by the user.", 54, 2,
	      stdout);
  PrintVADesc(1, "ntqs", 11, "FrameRate", 14, "The rate of charge density "
	      "sampling; the number of steps between successive snapshots to "
	      "determine the solvent reaction field potential (SRFP).  "
	      "Default 1000 (the time step set in the &cntrl or &ipolq "
	      "namelists should factor into the ntqs setting).", 54, 2,
	      stdout);
  PrintVADesc(1, "nqframe", 11, "FrameCount", 14, "The number of frames used "
	      "to compose the SRFP.  Default 10 (this is too low for most "
	      "solvent environments).", 54, 2, stdout);
  PrintVADesc(1, "nsteqlim", 11, "EqStepCount", 14, "The number of steps used "
	      "to equilibrate the system before charge density collection "
	      "begins.  Use this part of the simulation to buffer against any "
	      "artifacts that might arise from suddenly freezing the solute "
	      "in place.  Default 10000.", 54, 2, stdout);
  PrintVADesc(1, "nblock", 11, "Blocks", 14, "The number of blocks into which "
	      "the simulation shall be divided for the purpose of estimating "
	      "the convergence of the electrostatic potential.  Default 4.",
	      54, 2, stdout);
  PrintVADesc(1, "verbose", 11, "Verbose", 14, "Default 0; set to 1 to "
	      "activate step-by-step progress updates printed to the terminal "
	      "window.  Useful for monitoring short runs to ensure that the "
	      "input successfully completes the SRFP calculation.", 54, 2,
	      stdout);
  PrintVADesc(1, "econv", 11, "EConverge", 14, "Convergence tolerance for the "
	      "SRFP (convergence checking is not yet implemented).", 54, 2,
	      stdout);
  PrintVADesc(1, "nqshell", 11, "QShellCount", 14, "The number of additional "
	      "shells of charge to place around the system in order to "
	      "approximate the SRFP due to infinite electrostatics in the "
	      "confines of an isolated system.  Maximum (and default) is 3, "
	      "minimum is 1.", 54, 2, stdout);
  PrintVADesc(1, "nvshell", 11, "VShellCount", 14, "The number of shells "
	      "around each atom on which the exact SRFP due to infinite "
	      "electrostatics shall be calculated.  Maximum (and default) is "
	      "3, minimum is 1.", 54, 2, stdout);
  PrintVADesc(1, "nqphpt", 11, "QSpherePts", 14, "In order to generate the "
	      "surface charges that will help in approximating the SRFP, "
	      "this number of points is placed equidistant on a sphere.  The "
	      "sphere is then rotated randomly and expanded to the radii "
	      "indicated by qshell[1,2,3,x].  All points that are on the "
	      "sphere due to one atom but within the sphere projected by "
	      "another atom are deleted, until only points on the proper "
	      "surface remain.  Default 100.", 54, 2, stdout);
  PrintVADesc(1, "nvphpt", 11, "VSpherePts", 14, "Similar to nqphpt above, "
	      "but for the shells of SRFP evaluation points.  Default 20.", 54,
	      2, stdout);
  PrintVADesc(1, "qshell1", 11, "ExpQBoundary", 14, "The distance at which to "
	      "locate the first surface charges, and to stop collecting "
	      "charges explicitly from the simulation's non-solute (that is, "
	      "solvent) atoms.  Default 5.0.", 54, 2, stdout);
  PrintVADesc(1, "qshell[2,3]", 11, "QShell[2,3]", 14, "The distance at which "
	      "to locate the second and third shells of boundary charges.  If "
	      "engaged, each shell must be located successively further from "
	      "the solute than the previous one.", 54, 2, stdout);
  PrintVADesc(1, "qshellx", 11, "QShellX", 14, "The distance at which to "
	      "locate an interior shell of charges, and the minimum distance "
	      "at which to take charges explicitly from the simulation.  Use "
	      "this feature if there is a risk that QM basis functions might "
	      "latch onto the simulation's point charges in absence of any "
	      "exclusion effects, and thus distort the wavefunction.", 54, 2,
	      stdout);
  PrintVADesc(1, "vshell[1-3]", 11, "VShell[1-3]", 14, "The distances at "
	      "which to locate additional shells of exact SRFP evaluation "
	      "points.  The SRFP is always evaluated, exactly, at the solute "
	      "atom sites.", 54, 2, stdout);
  PrintVADesc(1, "dt", 11, "TimeStep", 14, "The simulation time step.  This "
	      "is read in the &ipolq namelist just as if it were present in "
	      "the &cntrl namelist, but a value specified in &ipolq overrides "
	      "the &cntrl setting.  Default is 0.001ps, set in &cntrl.", 54, 2,
	      stdout);
  PrintVADesc(1, "minqwt", 11, "MinQWeight", 14, "The stiffness of harmonic "
	      "restraint by which to restrain fitted shell charges to zero.  "
	      "Default 0.01.", 54, 2, stdout);
  PrintVADesc(1, "modq", 11, "ModifyQ", 14, "When IPolQ is applied, it is "
	      "appropriate to hyper-polarize certain molecules in the SRFP "
	      "calculation.  This variable may be specified as many times as "
	      "necessary, followed by an ambmask string and a real number "
	      "indicating the new charges to be assigned to all atoms in the "
	      "mask.  For example, fixed-charge water models should have "
	      "their dipoles increased by an amount equal to the original "
	      "model's dipole less 1.85 (the dipole of water in vacuum).", 54,
	      2, stdout);
  PrintVADesc(1, "prepqm", 11, "QuantumPrep", 14, "Preparatory call for QM "
	      "calculations.  This variable may be specified as many times as "
	      "necessary.  Each of these calls will be issued, in the order "
	      "specified, before executing quantum calculations.", 54, 2,
	      stdout);
  PrintVADesc(1, "postqm", 11, "QuantumClean", 14, "Post-processing calls for "
	      "QM calculations.  Similar to prepqm directives, called after "
	      "QM calculations have been completed.", 54, 2, stdout);
  PrintVADesc(1, "qmprog", 11, "QMPackage", 14, "The quantum package to use.  "
	      "Supported packages are \"gaussian\" and \"orca\".", 54, 2,
	      stdout);
  PrintVADesc(1, "qmpath", 11, "QMPath", 14, "Path to the primary QM "
	      "executable.  This path will be tested, taking into account "
	      "prepqm calls, to be sure that the executable exists prior to "
	      "running the SRFP calculation.", 54, 2, stdout);
  PrintVADesc(1, "qmcomm", 11, "QMInputFile", 14, "The base name of the QM "
	      "input file.  Vacuum and condensed-phase versions will be "
	      "written with extensions 'vacu' and 'solv', respectively.  "
	      "Default 'IPolQinp'.", 54, 2, stdout);
  PrintVADesc(1, "maxcore", 11, "MaxMemory", 14, "The maximum memory that can "
	      "be allocated to arrays for quantum calculations with Orca, or "
	      "the maximum total memory that can be allocated for "
	      "calculations with Gaussian.", 54, 2, stdout);
  PrintVADesc(1, "qmresult", 11, "QMOutputFile", 14, "The base name of the QM "
	      "output file, which is given similar extensions to the input "
	      "file.  Default 'IPolQout'.", 54, 2, stdout);
  PrintVADesc(1, "ptqfi", 11, "PointQFile", 14, "The name of the point "
	      "charges file referenced by orca for including the SRFP into "
	      "the condensed-phase calculation.", 54, 2, stdout);
  PrintVADesc(1, "qmflag", 11, "QMSignal", 14, "The name of the file used to "
	      "signal slave processes that the QM calculations launched by "
	      "the master are complete.  Default '.mdgx.finqm'.", 54, 2,
	      stdout);
  PrintVADesc(1, "qmlev", 11, "QMTheory", 14, "The level of QM theory to "
	      "use.  Default MP2.", 54, 2, stdout);
  PrintVADesc(1, "basis", 11, "QMBasis", 14, "The QM basis set to use.  "
	      "Default cc-pvTZ.", 54, 2, stdout);
  PrintVADesc(1, "scrdir", 11, "WorkDirectory", 14, "The scratch directory "
	      "to use during QM calculations.  Useful to reduce NFS load.  If "
	      "the directory exists, it will be used but not destroyed "
	      "following each QM calculation.  If the directory does not "
	      "exist at the start of the run, it will be created and later "
	      "destroyed.", 54, 2, stdout);
  PrintVADesc(1, "rqminp", 11, "KeepQMInput", 14, "Directive to retain QM "
	      "input files after the run.  Default 0 (OFF).", 54, 2, stdout);
  PrintVADesc(1, "rqmchk", 11, "KeepQMCheckPt", 14, "Directive to retain QM "
	      "checkpoint file(s) after the run.  Default 0 (OFF).", 54, 2,
	      stdout);
  PrintVADesc(1, "rqmout", 11, "KeepQMOutput", 14, "Directive to retain QM "
	      "output files after the run.  Default 0 (OFF).", 54, 2, stdout);
  PrintVADesc(1, "rcloud", 11, "KeepQCloud", 14, "Directive to retain the "
	      "solvent charge density cloud file after the run.  Default 0 "
	      "(OFF).", 54, 2, stdout);
  PrintVADesc(1, "checkex", 11, "CheckExist", 14, "Activates safety checks "
	      "for the existence of QM executables (including electrostatic "
	      "potential calculators) called at the start of the run.  These "
	      "checks attempt to take into account user-specified preparatory "
	      "directives (see prepqm above).  Default 1 (ON).  Set to zero "
	      "to disable this safeguard, for instance if the checks cannot "
	      "find the executables but the preparatory directives, when "
	      "fully implemented, are known to result in success.", 54, 2,
	      stdout);
  PrintVADesc(1, "unx", 11, "UElecXBin", 14, "The number of grid points on "
	      "which to evaluate the electrostatic potential, in the X "
	      "direction.  Grid dimensions in Y and Z are set by similar "
	      "variables.", 54, 2, stdout);
  PrintVADesc(1, "uhx", 11, "UElecXSpc", 14, "The grid spacing of the "
	      "electrostatic potential grid in the X direction.  The grid is "
	      "always rectilinear.  Spacings in Y and Z are set by similar "
	      "variables.", 54, 2, stdout);
  PrintVADesc(1, "cengrid", 11, "CenterGrid", 14, "Directive to center the "
	      "electrostatic potential grid on the location of the molecule "
	      "stored in mdgx.  The default behavior varies with each quantum "
	      "package: 'orca' activates centering on the molecule whereas "
	      "'gaussian' calls for centering on the origin, as Orca does not "
	      "reposition the molecule in its output but Gaussian will place "
	      "the molecule in a 'Standard Orientation' and leave it there "
	      "in the output and checkpoint files used for electrostatic "
	      "potential calculations.", 54, 2, stdout);
  PrintVADesc(1, "fmpath", 11, "FormChkPath", 14, "Path to the program called "
	      "for converting binary checkpoint files into formatted "
	      "checkpoint files.  Needed only if the QM program is 'gaussian'"
	      ".", 54, 2, stdout);
  PrintVADesc(1, "uvpath", 11, "UEvalPath", 14, "Path to the program called "
	      "for evaluating the electrostatic potential grid.", 54, 2,
	      stdout);
  PrintVADesc(1, "grid", 11, "GridFile", 14, "Base name of the electrostatic "
	      "potential grid to be written.  As with QM input and output, "
	      "this base name is appended 'vacu' or 'solv' for vacuum and "
	      "condensed-phase calculations.", 54, 2, stdout);
}

/***=======================================================================***/
/*** PrintAttributions: users are provided with a straightforward means of ***/
/***                    seeing which aspects of this code were taken from  ***/
/***                    other sources.  Attributions are also provided in  ***/
/***                    comments to the source code, where appropriate.    ***/
/***=======================================================================***/
void PrintAttributions()
{
  PrintSplash(stdout);
  PrintParagraph("Attributions:", 79, stdout);
  PrintParagraph("- Implementation of the SETTLE algorithm (J. Comput. Chem. "
		 "13:952-966, 1992) was adapted from the NAMD program source, "
		 "v2.6, developed by the Theoretical and Computational "
		 "Biophysics Group in the Beckman Institute for Advanced "
		 "Science and Technology at the University of Illinois at "
		 "Urbana-Champaign.", 79, stdout);
  PrintParagraph("- Implementation of the Smooth Particle Mesh Ewald "
		 "algorithm (J. Chem. Phys. 103, 8577-8593, 1995) includes "
		 "code developed by Thomas A. Darden for optimization of the "
		 "convolution kernel.", 79, stdout);
  PrintParagraph("- Implementation of the Langevin dynamics integration "
		 "algorithm (Biopolymers 32, 523-535, 1992) is adapted from "
		 "the sff program also distributed with AmberTools.", 79,
		 stdout);
  PrintParagraph("- Dr. Robert E. Duke is acknowledged for outstanding advice "
		 "and insights into the problem of efficient and scalable "
		 "molecular dynamics.  The mdgx program would not exist "
		 "without his support.", 79, stdout);
}
