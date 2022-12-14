#ifndef INC_DATAIO_EVECS_H
#define INC_DATAIO_EVECS_H
#include "DataIO.h"
/// Read write evecs (eigenmodes) data file.
class DataIO_Evecs : public DataIO {
  public:
    DataIO_Evecs();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Evecs(); }
    static void ReadHelp();
    int processReadArgs(ArgList &);
    int ReadData(std::string const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList &)                     { return 0; }
    int WriteData(std::string const&,DataSetList const&);
    int WriteData2D(std::string const&, DataSetList const&) { return 1; }
    int WriteData3D(std::string const&, DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&);
  private:
    static const char* MatrixOutputString(DataSet::scalarType);
    int ibeg_, iend_;
    bool hasIend_;
};
#endif
