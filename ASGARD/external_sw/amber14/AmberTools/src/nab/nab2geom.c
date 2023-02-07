#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#if AVS
#include <avs/avs.h>
#include <avs/geom.h>
#endif

#include "nab.h"

typedef	int	BOND_T[ 2 ];

#if !AVS	

char	*nab2geom( MOLECULE_T *mol )
{

	return( NULL );
}

#else

GEOMobj	*nab2geom( MOLECULE_T * );
static	void	setcolor( POINT_T, ATOM_T * );

GEOMobj	*nab2geom( mol )
{
	GEOMobj	*obj = NULL;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	int	r, tr;
	int	ri, rj;
	int	a, ta, b, tb, c;
	int	ai, aj, tp;
	float	xm, ym, zm;
	int	*a_res = NULL;
	int	*a_off = NULL;
	ATOM_T	**a_index = NULL;
	BOND_T	*bonds = NULL;
	POINT_T	*a_pos = NULL;
	POINT_T	*a_colors = NULL;
	ATOM_T	*api, *apj;
	EXTBOND_T	*ebp;

	for( tr = ta = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		tr += sp->s_nresidues;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			ta += res->r_natoms;
		}
	}
	tb = 1.5 * ta;

	if( ( a_res = ( int * )malloc( tr * sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "nab2geom: can't alloc a_res\n" );
		goto clean_up;
	}
	if( ( a_off = ( int * )malloc( tr * sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "nab2geom: can't alloc a_res\n" );
		goto clean_up;
	}
	if( ( a_index = ( ATOM_T ** )malloc( ta * sizeof( ATOM_T * ) ) )
		== NULL ){
		fprintf( stderr, "nab2geom: can't alloc a_index\n" );
		goto clean_up;
	}
	if( ( bonds = ( BOND_T * )malloc( tb * sizeof( BOND_T ) ) )
		== NULL ){
		fprintf( stderr, "nab2geom: can't alloc bonds\n" );
		goto clean_up;
	}

	for( tr = ta = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++, tr++ ){
			a_off[ tr ] = ta;
			res = sp->s_residues[ r ];
			a_res[ tr ] = ta;
			for( a = 0; a < res->r_natoms; a++, ta++ )
				a_index[ ta ] = &res->r_atoms[ a ];
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

	if( ( a_pos = ( POINT_T * )malloc( 4 * tb * sizeof( POINT_T ) ) )
		== NULL ){
		fprintf( stderr, "nab2geom: can't alloc a_pos\n" );
		goto clean_up;
	}

	if( ( a_colors = ( POINT_T * )malloc( 4 * tb * sizeof( POINT_T ) ) )
		== NULL ){
		fprintf( stderr, "nab2geom: can't alloc a_colors\n" );
		goto clean_up;
	}

	for( tp = 0, b = 0; b < tb; b++, tp += 4 ){
		ai = bonds[ b ][ 0 ];
		aj = bonds[ b ][ 1 ];
		api = a_index[ ai ];
		apj = a_index[ aj ];
		xm = 0.5 * ( api->a_pos[ 0 ] + apj->a_pos[ 0 ] );
		ym = 0.5 * ( api->a_pos[ 1 ] + apj->a_pos[ 1 ] );
		zm = 0.5 * ( api->a_pos[ 2 ] + apj->a_pos[ 2 ] );
		a_pos[ tp ][ 0 ] = api->a_pos[ 0 ];
		a_pos[ tp ][ 1 ] = api->a_pos[ 1 ];
		a_pos[ tp ][ 2 ] = api->a_pos[ 2 ];
		a_pos[ tp + 1 ][ 0 ] = xm;
		a_pos[ tp + 1 ][ 1 ] = ym;
		a_pos[ tp + 1 ][ 2 ] = zm;
		a_pos[ tp + 2 ][ 0 ] = xm;
		a_pos[ tp + 2 ][ 1 ] = ym;
		a_pos[ tp + 2 ][ 2 ] = zm;
		a_pos[ tp + 3 ][ 0 ] = apj->a_pos[ 0 ];
		a_pos[ tp + 3 ][ 1 ] = apj->a_pos[ 1 ];
		a_pos[ tp + 3 ][ 2 ] = apj->a_pos[ 2 ];
		setcolor( a_colors[ tp ], api );
		setcolor( a_colors[ tp + 1 ], api );
		setcolor( a_colors[ tp + 2 ], apj );
		setcolor( a_colors[ tp + 3 ], apj );
	}

	if( !( obj = GEOMcreate_obj( GEOM_POLYTRI, GEOM_NULL ) ) ){
		fprintf( stderr, "nab2geom: can't create obj\n" );
		goto clean_up;
	}
	GEOMadd_disjoint_line( obj,
		a_pos, a_colors, tp, GEOM_DONT_COPY_DATA );

clean_up : ;
	if( bonds )
		free( bonds );
	if( a_index )
		free( a_index );
	if( a_off )
		free( a_off );
	if( a_res )
		free( a_res );

	return( obj );
}

static	void	setcolor( POINT_T cval, ATOM_T *ap )
{
	int	atype;

	atype = *ap->a_name;
	if( isdigit( atype ) )
		atype = ap->a_name[ 1 ];
	switch( atype ){
	case 'c' :
	case 'C' :
		cval[ 0 ] = 0.0;
		cval[ 1 ] = 1.0;
		cval[ 2 ] = 0.0;
		break;
	case 'o' :
	case 'O' :
		cval[ 0 ] = 1.0;
		cval[ 1 ] = 0.0;
		cval[ 2 ] = 0.0;
		break;
	case 'n' :
	case 'N' :
		cval[ 0 ] = 0.0;
		cval[ 1 ] = 1.0;
		cval[ 2 ] = 1.0;
		break;
	case 'h' :
	case 'H' :
		cval[ 0 ] = 1.0;
		cval[ 1 ] = 1.0;
		cval[ 2 ] = 1.0;
		break;
	case 's' :
	case 'S' :
		cval[ 0 ] = 1.0;
		cval[ 1 ] = 1.0;
		cval[ 2 ] = 0.0;
		break;
	case 'p' :
	case 'P' :
		cval[ 0 ] = 1.0;
		cval[ 1 ] = 1.0;
		cval[ 2 ] = 0.0;
		break;
	default :
		cval[ 0 ] = 1.0;
		cval[ 1 ] = 0.0;
		cval[ 2 ] = 1.0;
		break;
	}
}

#endif	/*  !AVS  */
