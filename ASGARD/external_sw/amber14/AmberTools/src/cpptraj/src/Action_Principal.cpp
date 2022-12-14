#include "Action_Principal.h"
#include "CpptrajStdio.h"
#include "Matrix_3x3.h"

// CONSTRUCTOR
Action_Principal::Action_Principal() :
  doRotation_(false),
  useMass_(false),
  debug_(0),
  outfile_(0),
  vecData_(0),
  valData_(0)
{ }

void Action_Principal::Help() {
  mprintf("\t[<mask>] [dorotation] [out <filename>] [name <dsname>]\n"
          "  Calculate principal axes of atoms in <mask>. Align the system along\n"
          "  principal axes if 'dorotation' specified.\n");
}

// Action_Principal::Init()
Action::RetType Action_Principal::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Keywords
  std::string dsname = actionArgs.GetStringKey("name");
  doRotation_ = actionArgs.hasKey("dorotation");
  // CPPTRAJ always uses mass no matter what this keyword says.
  useMass_ = actionArgs.hasKey("mass");
  std::string filename = actionArgs.GetStringKey("out");
  if (!doRotation_ && filename.empty() && dsname.empty()) {
    mprinterr("Error: At least one of 'dorotation', 'out <filename>', or 'name <dsname>' must be specified.\n");
    return Action::ERR;
  }
  // Masks
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  // Set up data
  if (!dsname.empty()) {
     vecData_ = (DataSet_Mat3x3*)DSL->AddSetAspect(DataSet::MAT3X3, dsname, "evec");
     valData_ = (DataSet_Vector*)DSL->AddSetAspect(DataSet::VECTOR, dsname, "eval");
     if (vecData_ == 0 || valData_ == 0) return Action::ERR;
  }

  mprintf("    PRINCIPAL:");
  if (!filename.empty()) {
    outfile_ = DFL->AddCpptrajFile(filename, "Eigenvectors/Eigenvalues");
    if (outfile_ == 0) return Action::ERR;
    mprintf(" output eigenvectors/eigenvalues to %s,", outfile_->Filename().full());
  }
  if (doRotation_)
    mprintf(" with rotation by");
  else
    mprintf(" without rotation by");
  if (useMass_)
    mprintf(" center of mass");
  else
    mprintf(" center of geometry");
  mprintf(", atoms selected by [%s]\n", mask_.MaskString());
  if (vecData_ != 0)
    mprintf("\tSaving eigenvectors to '%s' (in rows of 3x3 matrices).\n"
            "\tSaving eigenvalues to '%s'\n", vecData_->legend(), valData_->legend());

  return Action::OK;
}

// Action_Principal::Setup()
Action::RetType Action_Principal::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask(mask_)) return Action::ERR;

  if (mask_.None()) {
    mprintf("Warning: No atoms selected for %s [%s].\n",currentParm->c_str(), mask_.MaskString());
    return Action::ERR;
  }

  mprintf("\tSelected %i atoms.\n", mask_.Nselected());
  return Action::OK;
}

// Action_Principal::DoAction()
Action::RetType Action_Principal::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 Inertia;
  Vec3 Eval;

  currentFrame->CalculateInertia( mask_, Inertia );

  // NOTE: Diagonalize_Sort_Chirality places sorted eigenvectors in rows.
  Inertia.Diagonalize_Sort_Chirality( Eval, debug_ );
  if (outfile_ != 0) {
    int fn = frameNum+1; 
    outfile_->Printf("%i EIGENVALUES: %f %f %f\n%i EIGENVECTOR 0: %f %f %f\n%i EIGENVECTOR 1: %f %f %f\n%i EIGENVECTOR 2: %f %f %f\n", 
      fn, Eval[0], Eval[1], Eval[2],
      fn, Inertia[0], Inertia[1], Inertia[2],
      fn, Inertia[3], Inertia[4], Inertia[5],
      fn, Inertia[6], Inertia[7], Inertia[8]);
    //Eval.Print("PRINCIPAL EIGENVALUES");
    //Inertia.Print("PRINCIPAL EIGENVECTORS (Rows)");
  }
  if (vecData_ != 0) {
    vecData_->AddMat3x3( Inertia );
    valData_->AddVxyz( Eval );
  }
  
  // Rotate - since Evec is already transposed (eigenvectors
  // are returned in rows) just do plain rotation to affect an
  // inverse rotation.
  if (doRotation_)
    currentFrame->Rotate( Inertia );

  return Action::OK;
}
