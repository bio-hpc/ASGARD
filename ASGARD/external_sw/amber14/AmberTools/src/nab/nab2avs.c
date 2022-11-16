#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#if AVS 
#include "avs.h"
#include "field.h"
#endif

#include "nab.h"

typedef	int	BOND_T[ 2 ];

#if !AVS

int	nab2avs( MOLECULE_T *mol )
{

	return( 0 );
}

#else

int	nab2avs( MOLECULE_T *, AVSfield_char **, AVSfield_int **, AVSfield_float ** );
MOLECULE_T	*avs2nab( AVSfield_char *, AVSfield_int *, AVSfield_float * );

int	nab2avs( MOLECULE_T *mol,
	AVSfield_char **AtomNames, AVSfield_int **AtomBonds, AVSfield_float **AtomXYZP )
{
	STRAND_T	*sp;
	RESIDUE_T	*res;
	int	r, tr;
	int	ri, rj;
	int	a, ta, b, tb, c;
	int	ai, aj, tp;
	int	dims0[ 1 ];
	int	*a_res = NULL;
	int	*a_off = NULL;
	int	*a_rnum = NULL;
	ATOM_T	**a_index = NULL;
	BOND_T	*bonds = NULL;
	POINT_T	*a_pos = NULL;
	ATOM_T	*ap, *api, *apj;
	EXTBOND_T	*ebp;
	unsigned char	 *np;
	int	*bp;
	float	*xp;

	for( tr = ta = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		tr += sp->s_nresidues;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			ta += res->r_natoms;
		}
	}
	tb = 1.5 * ta;

	if( ( a_res = ( int * )malloc( tr * sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "nab2avs: can't alloc a_res\n" );
		goto clean_up;
	}
	if( ( a_off = ( int * )malloc( tr * sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "nab2avs: can't alloc a_res\n" );
		goto clean_up;
	}
	if( ( a_index = ( ATOM_T ** )malloc( ta * sizeof( ATOM_T * ) ) )
		== NULL ){
		fprintf( stderr, "nab2avs: can't alloc a_index\n" );
		goto clean_up;
	}
	if( ( a_rnum = ( int * )malloc( ta * sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "nab2avs: can't alloc a_rnum\n" );
		goto clean_up;
	}
	if( ( a_pos = ( POINT_T * )malloc( ta * sizeof( POINT_T ) ) ) == NULL ){
		fprintf( stderr, "nab2avs: can't alloc a_pos\n" );
		goto clean_up;
	}
	if( ( bonds = ( BOND_T * )malloc( tb * sizeof( BOND_T ) ) )
		== NULL ){
		fprintf( stderr, "nab2avs: can't alloc bonds\n" );
		goto clean_up;
	}

	for( tr = ta = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++, tr++ ){
			a_off[ tr ] = ta;
			res = sp->s_residues[ r ];
			a_res[ tr ] = ta;
			for( a = 0; a < res->r_natoms; a++, ta++ ){
				a_index[ ta ] = &res->r_atoms[ a ];
				a_rnum[ ta ] = tr;
			}
		}
	}

	for( tr = ta = tb = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			for( ai = 0; ai < res->r_natoms; ai++ ){
				api = &res->r_atoms[ ai ];
				for( c = 0; c < api->a_nconnect; c++ ){
					aj = api->a_connect[ c ];
					if( ai < aj ){
						bonds[ tb ][ 0 ] = ai + ta;
						bonds[ tb ][ 1 ] = aj + ta;
						tb++;
					}
				}
			}
			for( ebp = res->r_extbonds; ebp; ebp = ebp->eb_next ){
				if( ( rj = ebp->eb_rnum ) < r + 1 )
					continue;
				ai = a_off[ tr + r ];
				aj = a_off[ tr + rj - 1 ];
				bonds[ tb ][ 0 ] = ebp->eb_anum + ai -1;
				bonds[ tb ][ 1 ] = ebp->eb_ranum + aj -1;
				tb++;
			} 
			ta += res->r_natoms;
		}
		tr += sp->s_nresidues;
	}

	if( *AtomNames )
		AVSfield_free( *AtomNames );
	dims0[ 0 ] = ta * 16;
	*AtomNames = ( AVSfield_char * )AVSdata_alloc(
		"field 1D 1-space 1-vector uniform byte", dims0 );

	if( *AtomBonds )
		AVSfield_free( *AtomBonds );
	dims0[ 0 ] = tb;
	*AtomBonds = ( AVSfield_int * )AVSdata_alloc(
		"field 1D 1-space 2-vector uniform integer", dims0 );

	if( *AtomXYZP )
		AVSfield_free( *AtomXYZP );
	dims0[ 0 ] = ta;
	*AtomXYZP = ( AVSfield_float * )AVSdata_alloc(
		"field 1D 1-space 4-vector uniform float", dims0 );
/*
	AVSfield_set_labels( *AtomXYZ, "X;Y;Z;float1", ";" );
*/
	np = ( *AtomNames )->data;
	for( a = 0; a < ta; a++ ){
		ap = a_index[ a ];
		strcpy( np, ap->a_atomname );
		strcpy( &np[ 5 ], ap->a_residue->r_resname );
		sprintf( &np[ 10 ], "%d",  a_rnum[ a ] );
		np += 16;
	}
	bp = ( *AtomBonds )->data;
	for( b = 0; b < tb; b++ ){
		bp[ 0 ] = bonds[ b ][ 0 ];
		bp[ 1 ] = bonds[ b ][ 1 ];
		bp += 2;
	}
	xp = ( *AtomXYZP )->data;
	for( a = 0; a < ta; a++ ){
		ap = a_index[ a ];
		xp[ 0 ] = ap->a_pos[ 0 ];
		xp[ 1 ] = ap->a_pos[ 1 ];
		xp[ 2 ] = ap->a_pos[ 2 ];
		xp[ 3 ] = ap->a_float1;
		/*xp[ 4 ] = ap->a_float2;*/
		xp += 4;
	}

clean_up : ;
	if( bonds )
		free( bonds );
	if( a_rnum )
		free( a_rnum );
	if( a_index )
		free( a_index );
	if( a_off )
		free( a_off );
	if( a_res )
		free( a_res );

	return( 1 );
}

MOLECULE_T	*avs2nab( AVSfield_char *AtomNames, AVSfield_int *AtomBonds,
	AVSfield_float *AtomXYZP )
{
	MOLECULE_T	*mp;
	STRAND_T	*sp;
	RESIDUE_T	res;
	RESIDUE_T	*rpi, *rpj;
	int	r, lr, nres;
	int	ri, rj;
	int	a, natoms, acnt, tacnt; 
	int	ai, aj;
	int	b, nbonds;
	int	*a_res;
	int	*r_atom1;
	BOND_T	bonds[ 150 ];
	ATOM_T	atoms[ 100 ];
	ATOM_T	*ap, *api, *apj;
	EXTBOND_T	*ebp;
	unsigned char	 *np;
	int	*bp;
	float	*xp;
	
	natoms = MAXX( AtomXYZP );
	nbonds = MAXX( AtomBonds );

	bzero( &res, sizeof( res ) );
	res.r_atoms = atoms;
	mp = ( MOLECULE_T * )newmolecule();
	addstrand( mp, "str_1" );

	np = AtomNames->data;
	nres = atoi( &np[ ( natoms - 1 ) * 16 + 10 ] );

	if( !( a_res = ( int * )malloc( natoms * sizeof( int ) ) ) ){
		fprintf( stderr, "avs2nab: can't allocate a_res\n" );
		return( NULL );
	}
	if( !( r_atom1 = ( int * )malloc( nres * sizeof( int ) ) ) ){
		fprintf( stderr, "avs2nab: can't allocate r_atom1\n" );
		return( NULL );
	}

	xp = AtomXYZP->data;
	r_atom1[ 0 ] = tacnt = 0;
	for( lr = -1, acnt = 0, a = 0; a < natoms; a++ ){
		if( ( r = atoi( &np[ 10 ] ) ) != lr ){
			if( acnt > 0 ){
				res.r_natoms = acnt;
				addresidue( mp, "str_1", &res );
				tacnt += acnt;
				r_atom1[ r ] = tacnt;
			}
			strcpy( res.r_resname, &np[ 5 ] );
			acnt = 0;
		}
		ap = &atoms[ acnt ];
		bzero( ap, sizeof( ATOM_T ) );
		strcpy( ap->a_atomname, np );
		ap->a_residue = &res;
		ap->a_pos[ 0 ] = xp[ 0 ];
		ap->a_pos[ 1 ] = xp[ 1 ];
		ap->a_pos[ 2 ] = xp[ 2 ];
		if( AtomXYZP->veclen > 3 )
			ap->a_float1 = xp[ 3 ];
		xp += AtomXYZP->veclen;
		acnt++;
		a_res[ a ] = r;
		lr = r;
		np += 16;
	}
	if( acnt > 0 ){
		res.r_natoms = acnt;
		addresidue( mp, "str_1", &res );
	}

	sp = mp->m_strands;
	bp = AtomBonds->data;
	for( b = 0; b < nbonds; b++ ){
		ai = bp[ 0 ];
		ri = a_res[ ai ];
		aj = bp[ 1 ];
		rj = a_res[ aj ];
		if( ri == rj ){
			rpi = sp->s_residues[ ri ];
			ai -= r_atom1[ ri ];
			aj -= r_atom1[ ri ];
			api = &rpi->r_atoms[ ai ];
			apj = &rpi->r_atoms[ aj ];
			api->a_connect[ api->a_nconnect ] = aj;
			api->a_nconnect++;
			apj->a_connect[ apj->a_nconnect ] = ai;
			apj->a_nconnect++;
			apj = &rpi->r_atoms[ aj ];
		}else{
			rpi = sp->s_residues[ ri ];
			rpj = sp->s_residues[ rj ];
			ai -= r_atom1[ ri ];
			aj -= r_atom1[ rj ];
			ebp = ( EXTBOND_T * )malloc( sizeof( EXTBOND_T ) );
			ebp->eb_anum = ai + 1;
			ebp->eb_rnum = rj + 1;
			ebp->eb_ranum = aj + 1;
			ebp->eb_next = rpi->r_extbonds;
			rpi->r_extbonds = ebp;
			ebp = ( EXTBOND_T * )malloc( sizeof( EXTBOND_T ) );
			ebp->eb_anum = aj + 1;
			ebp->eb_rnum = ri + 1;
			ebp->eb_ranum = ai + 1;
			ebp->eb_next = rpj->r_extbonds;
			rpj->r_extbonds = ebp;
		}
		bp += 2;
	}
	
	free( r_atom1 );
	free( a_res );

	return( mp );
}

#endif
