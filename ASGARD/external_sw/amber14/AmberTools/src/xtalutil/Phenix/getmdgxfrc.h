#ifndef _thisheaderfile_H_
#define _thisheaderfile_H_


extern "C"
{
    #include "mdgx.h"

    
    int getmdgxfrc(const double *PhenixCoords, double* target,
                   double * gradients, uform* U, trajcon* tjptr, mdsys* MDptr);
    
    trajcon CreateTrajCon();
    uform LoadTopology(const char *tpname, trajcon* tj);
    mdsys CreateMDSys(const char *crdname, uform* U);
    void DestroyMDSys(mdsys* mdsysptr);
    void DestroyUform(uform* uptr, mdsys* mdsysptr);
    
}


#endif

