#ifndef INC_ACTION_ROTATE_H
#define INC_ACTION_ROTATE_H
#include "Action.h"
#include "DataSet_Mat3x3.h"
class Action_Rotate : public Action {
  public:
    Action_Rotate();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Rotate(); }
    static void Help();
  private:
    Matrix_3x3 RotMatrix_;
    AtomMask mask_;
    DataSet_Mat3x3* rmatrices_;
    bool inverse_;

    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}
};
#endif
