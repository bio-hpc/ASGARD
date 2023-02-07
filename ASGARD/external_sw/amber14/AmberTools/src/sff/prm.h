#ifndef  PRM_H
#define  PRM_H

/*  All these declarations should be moved to sff.h after sff.h and nab*.h */
/*  are cleaned up to eliminate duplication that causes compilation failures */
/*  ../sff/sff.h:25: error: redefinition of typedef ‘INT_T’ */

/*  prmtop routine public interfaces: */
INT_T    free_prm( PARMSTRUCT_T* prm );
PARMSTRUCT_T*  rdparm( STRING_T* name );
INT_T    wrparm( PARMSTRUCT_T* prm, STRING_T* name );

/*  mme routine public interfaces: */
INT_T    mme_init_sff( PARMSTRUCT_T*, INT_T*, INT_T*, REAL_T*, FILE_T* );

#endif
