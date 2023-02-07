#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"
#include "errormsg.h"
#include "geometry.h"
#include "machine.h"
#include "y.tab.h"

#define	MAXGEOMS	20
int	n_geoms = 0;
GEOM_T	geoms[ MAXGEOMS ];
GEOM_T	*cur_gp = NULL;
STR_T	*cur_sp = NULL;
RES_T	*cur_rp = NULL;

#define	R_ATTR_GEOM_T	0
#define	R_ATTR_STR	1
#define	R_ATTR_RES	2

typedef	struct	r_attr	{
	char	*in;
	char	*out;
	char	*reslib;
} R_ATTR;

int	r_attr_lev;
R_ATTR	r_attr[ 3 ];

int	r_list_lev;
RES_T	*r_list[ 2 ];
int	n_r_list_lev;
int	n_r_list[ 2 ];

void	dump_geom( FILE *fp )
{
	int	i, j;
	STR_T	*sp;
	RES_T	*rp;

	fprintf( fp, "Geometry:\t%s,", cur_gp->g_name );
	fprintf( fp, "  %2d strands\n", cur_gp->g_n_str );
	fprintf( fp, "  Helix :    " );
	for( i = 0; i < 6; i++ )
		fprintf( fp, "%8.3f", cur_gp->g_helix[ i ] );
	putc( '\n', fp );
	for( i = 0; i < cur_gp->g_n_str; i++ ){
		sp = &cur_gp->g_strs[ i ];
		fprintf( fp, "  Strand:\t%s, %2d residues",
			sp->s_name, sp->s_n_res );   
		if( sp->s_ref )
			fprintf( fp, ", reference" );
		putc( '\n', fp );
		fprintf( fp, "    Offset:  " );
		if( sp->s_offset ){
			for( j = 0; j < 6; j++ )
				fprintf( fp, "%8.3f", ( *sp->s_offset)[ j ] );
		}else
			fprintf( fp, "   none" );
		putc( '\n', fp );
		fprintf( fp, "    Orient:  " );
		if( sp->s_orient == G_PARALLEL )
			fprintf( fp, "   parallel" );
		else
			fprintf( fp, "   antiparallel" );
		putc( '\n', fp );
		fprintf( fp, "    Def Res: " );
		if( sp->s_defres == G_NONE )
			fprintf( fp, "   none" );
		else if( sp->s_defres == G_WATSONCRICK )
			fprintf( fp, "   watson/crick" );
		else
			fprintf( fp, "   copy" );
		putc( '\n', fp );
		for( j = 0; j < sp->s_n_res; j++ ){
			rp = &sp->s_res[ j ];
			fprintf( fp, "    Residue:\t%s\n", rp->r_reschar );
			fprintf( fp, "      Resname:  %s\n", rp->r_resname );
			fprintf( fp, "      In:       %s\n", rp->r_in );
			fprintf( fp, "      Out:      %s\n", rp->r_out );
			fprintf( fp, "      Reslib:   %s\n", rp->r_reslib );
		}
	}
}

void	mk_geom( int op, int val, NODE_T *left, NODE_T *right )
{
	R_ATTR	*rap;

	switch( op ){
	case GEOMETRY :
		if( val == G_START ){
			r_attr_lev = R_ATTR_GEOM_T;
			rap = &r_attr[ r_attr_lev ];
			rap->in = NULL;
			rap->out = NULL;
			rap->reslib = NULL;
			r_list_lev = R_ATTR_GEOM_T;
			r_list[ r_list_lev ] = NULL;
			n_r_list_lev = R_ATTR_GEOM_T;
			n_r_list[ n_r_list_lev ] = 0;
			new_geometry( left, right );
		}else{
			complete_geometry();
		}
		break;
	case HELIX :
		list2tform( cur_gp->g_helix, left ); 
		break;
	case RESLIB :
		if( r_attr[ r_attr_lev ].reslib )
			errormsg( TRUE, E_DUP_RESLIB );
		else
			r_attr[ r_attr_lev ].reslib =
				left->n_val.v_value.v_cval;
		break;
	case IN :
		if( r_attr[ r_attr_lev ].in )
			errormsg( TRUE, E_DUP_IN );
		else
			r_attr[ r_attr_lev ].in = left->n_val.v_value.v_cval;
		break;
	case OUT :
		if( r_attr[ r_attr_lev ].out )
			errormsg( TRUE, E_DUP_OUT );
		else
			r_attr[ r_attr_lev ].out = left->n_val.v_value.v_cval;
		break;
	case RESIDUE :
		if( val == G_START ){
			r_attr_lev = R_ATTR_RES;
			rap = &r_attr[ r_attr_lev ];
			rap->in = NULL;
			rap->out = NULL;
			rap->reslib = NULL;
			new_residue();
		}else{
			complete_residue();
			r_attr_lev = ( cur_sp != NULL ) ? 
				R_ATTR_STR : R_ATTR_GEOM_T;
			cur_rp = NULL;
		}
		break;
	case RESNAME :
		cur_rp->r_resname = left->n_val.v_value.v_cval;
		break;
	case RESCHAR :
		cur_rp->r_reschar = left->n_val.v_value.v_cval;
		break;
	case STRAND :
		if( val == G_START ){
			r_attr_lev = R_ATTR_STR;
			rap = &r_attr[ r_attr_lev ];
			rap->in = NULL;
			rap->out = NULL;
			rap->reslib = NULL;
			r_list_lev = R_ATTR_STR;
			r_list[ r_list_lev ] = NULL;
			n_r_list_lev = R_ATTR_STR;
			n_r_list[ n_r_list_lev ] = 0;
			new_strand( val, left, right );
		}else{
			complete_strand();
			r_attr_lev = R_ATTR_GEOM_T;
			r_list_lev = R_ATTR_GEOM_T;
			n_r_list_lev = R_ATTR_GEOM_T;
			cur_sp = NULL;
		}
		break;
	case OFFSET :
		if( left ){
			if( ( cur_sp->s_offset = 
				( TRANSFORM_T * )
				malloc( sizeof( TRANSFORM_T ) ) ) == NULL ){
				errormsg( TRUE, E_NOMEM_FOR_S, "offset" );
			}
			list2tform( cur_sp->s_offset, left );
		}
		break;
	case ORIENT :
		cur_sp->s_orient = val;
		break;
	case DEFRES :
		cur_sp->s_defres = val;
		break;
	default :
		errormsg_d( TRUE, E_UNEXPECTED_SYM_D, op );
		break;
	}
}

void	new_geometry( NODE_T *g_name, NODE_T *s_list )
{
	NODE_T	*lp;
	int	s;
	STR_T	*sp;
	SYMREC_T	*symp;

	if( findsym( g_name->n_val.v_value.v_cval ) )
		errormsg_s( TRUE, E_SYM_REDEF_S,
			g_name->n_val.v_value.v_cval );
	else if( ( symp = entersym( S_GLOBAL, g_name->n_val.v_value.v_cval,
		T_MOLECULE, C_GEOM, 0 ) ) == NULL )
		errormsg( TRUE, E_INTERNAL_ERROR );

	cur_gp = &geoms[ n_geoms++ ];
	cur_gp->g_name = g_name->n_val.v_value.v_cval;
	cur_gp->g_n_str = 0;
	for( lp = s_list; ; ){
		if( lp->n_sym == LIST ){
			if( lp->n_left->n_sym == IDENT )
				cur_gp->g_n_str++;
			lp = lp->n_right;
		}else if( lp->n_sym == IDENT ){
			cur_gp->g_n_str++;
			break;
		}
	}
	if( ( cur_gp->g_strs =
		( STR_T * )malloc( cur_gp->g_n_str * sizeof( STR_T ) ) ) ==
			NULL ){
		errormsg_s( TRUE, E_NOMEM_FOR_S, "strands" );
	}
	for( s = 0, sp = cur_gp->g_strs; s < cur_gp->g_n_str; s++, sp++ ){
		sp->s_ref = FALSE;
		sp->s_offset = NULL;
		sp->s_orient = G_UNDEF;
		sp->s_defres = G_UNDEF;
		sp->s_n_res = 0;
		sp->s_res = NULL;
	}
	for( lp = s_list, s = 0; ; ){
		if( lp->n_sym == LIST ){
			if( lp->n_left->n_sym == IDENT ){
				cur_gp->g_strs[ s ].s_name =
					lp->n_left->n_val.v_value.v_cval;
				s++;
			}
			lp = lp->n_right;
		}else if( lp->n_sym == IDENT ){
			cur_gp->g_strs[ s ].s_name = lp->n_val.v_value.v_cval;
			break;
		}
	}
}

void	complete_geometry( void )
{
	int	s;
	int	ref;
	STR_T	*sp;

	for( s = 0, ref = FALSE; s < cur_gp->g_n_str; s++ ){
		sp = &cur_gp->g_strs[ s ];
		if( sp->s_ref ){
			if( !ref )
				ref = TRUE;
			else{
				errormsg_s( TRUE, E_DUP_REF_STRAND_S,
					cur_gp->g_name );
			}
		}
		if( sp->s_n_res == 0 ){
			errormsg_2s( TRUE, E_STRAND_NOT_DEF_2S,
				sp->s_name, cur_gp->g_name );
		}
	}
	if( !ref )
		errormsg_2s( TRUE, E_GEOM_NO_REF_S, cur_gp->g_name );
}

void	new_strand( int val, NODE_T *id, int ref )
{
	int	s;
	STR_T	*sp;

	for( s = 0; s < cur_gp->g_n_str; s++, sp++ ){
		sp = &cur_gp->g_strs[ s ];
		if( strcmp( id->n_val.v_value.v_cval, sp->s_name ) == 0 ){
			cur_sp = sp;
			cur_sp->s_ref = ref == G_REFERENCE;
			cur_sp->s_offset = NULL;
			cur_sp->s_orient = G_UNDEF;
			cur_sp->s_defres = G_NONE;
			cur_sp->s_n_res = 0;
			return;
		}
	}
	errormsg_s( TRUE, E_MISSING_STRAND_S, id->n_val.v_value.v_cval );
}

void	complete_strand( void )
{
	int	nr;
	int	rs;
	RES_T	*rp, *rps, *rpg;
	char	msg[ 80 ];

	if( cur_sp->s_orient == G_UNDEF )
		errormsg_2s( TRUE, E_MISSING_STR_ITEM_2S,
			cur_sp->s_name, "orient" );

	for( nr = 0, rps = r_list[ r_list_lev ]; rps; rps = rps->r_next )
		nr++;
	for( rpg = r_list[ R_ATTR_GEOM_T ]; rpg; rpg = rpg->r_next ){
		for( rs = FALSE, rps = r_list[ r_list_lev ];
			rps; rps = rps->r_next ){
			if( strcmp( rps->r_reschar, rpg->r_reschar ) == 0 ){
				rs = TRUE;
				break;
			}
		}
		if( !rs )
			nr++;
	}
	if( ( rp = ( RES_T * )malloc( nr * sizeof( RES_T ) ) ) == NULL ){
		sprintf( msg, "residues of strand %s", cur_sp->s_name );
		errormsg_s( TRUE, E_NOMEM_FOR_S, msg );
	}
	cur_sp->s_n_res = nr;
	cur_sp->s_res = rp;

	for( nr = 0, rps = r_list[ r_list_lev ]; rps; rps = rps->r_next ){
		rp[ nr ].r_in = rps->r_in;
		rp[ nr ].r_out = rps->r_out;
		rp[ nr ].r_reslib = rps->r_reslib;
		rp[ nr ].r_reschar = rps->r_reschar;
		rp[ nr ].r_resname = rps->r_resname;
		nr++;
	}
	for( rpg = r_list[ R_ATTR_GEOM_T ]; rpg; rpg = rpg->r_next ){
		for( rs = FALSE, rps = r_list[ r_list_lev ];
			rps; rps = rps->r_next ){
			if( strcmp( rps->r_reschar, rpg->r_reschar ) == 0 ){
				rs = TRUE;
				break;
			}
		}
		if( !rs ){
			rp[ nr ].r_in = rpg->r_in;
			rp[ nr ].r_out = rpg->r_out;
			rp[ nr ].r_reslib = rpg->r_reslib;
			rp[ nr ].r_reschar = rpg->r_reschar;
			rp[ nr ].r_resname = rpg->r_resname;
			nr++;
		}
	}
}

void	new_residue( void )
{
	RES_T	*rp;

	if( ( rp = ( RES_T * )malloc( sizeof( RES_T ) ) ) == NULL )
		errormsg_s( TRUE, E_NOMEM_FOR_S, "residue" );
	rp->r_next = NULL;
	rp->r_in = NULL;
	rp->r_out = NULL;
	rp->r_reslib = NULL;
	rp->r_resname = NULL;
	rp->r_reschar = NULL;
	cur_rp = rp; 
}

void	complete_residue( void )
{
	int	cv, l;
	R_ATTR	*rap;
	RES_T	*rp, *rpl;

	for( l = R_ATTR_RES; l >= R_ATTR_GEOM_T; ){
		rap = &r_attr[ l ];
		if( rap->in ){
			if( cur_rp->r_in == NULL )
				cur_rp->r_in = rap->in;
		}
		if( rap->out ){
			if( cur_rp->r_out == NULL )
				cur_rp->r_out = rap->out;
		}
		if( rap->reslib ){
			if( cur_rp->r_reslib == NULL )
				cur_rp->r_reslib = rap->reslib;
		}
		if( l == R_ATTR_RES ){
			l = ( cur_sp ) ? l - 1 : l - 2;
		
		}else
			l--;
	}
	if( !cur_rp->r_reschar )
		errormsg( TRUE, E_RES_NOCHAR );
	if( !cur_rp->r_in )
		errormsg_2s( TRUE, E_MISSING_RES_ITEM_2S,
			cur_rp->r_reschar, "in" );
	if( !cur_rp->r_out )
		errormsg_2s( TRUE, E_MISSING_RES_ITEM_2S,
			cur_rp->r_reschar, "out" );
	if( !cur_rp->r_reslib )
		errormsg_2s( TRUE, E_MISSING_RES_ITEM_2S,
			cur_rp->r_reschar, "reslib" );
	if( !cur_rp->r_reschar )
		errormsg_2s( TRUE, E_MISSING_RES_ITEM_2S,
			cur_rp->r_reschar, "resname" );

	for( rpl = NULL, rp = r_list[ r_list_lev ]; rp;
			rpl = rp, rp = rp->r_next ){
		if( ( cv = strcmp( rp->r_reschar, cur_rp->r_reschar ) ) == 0 ){
			errormsg_s( TRUE, E_DUP_RES_S, cur_rp->r_reschar ); 
		}else if( cv > 0 )
			break;
	}
	cur_rp->r_next = rp;
	if( rpl )
		rpl->r_next = cur_rp;
	else
		r_list[ r_list_lev ] = cur_rp;
	n_r_list[ n_r_list_lev ]++;
}

void	list2tform( TRANSFORM tform, NODE_T *list )
{
	NODE_T	*lp;
	int	p;

	for( lp = list, p = 0; ; ){
		if( lp->n_sym == LIST ){
			if( lp->n_left->n_sym == INT_LIT ){
				tform[ p ] =
					lp->n_left->n_val.v_value.v_ival;
				p++;
			}else if( lp->n_left->n_sym == FLOAT_LIT ){
				tform[ p ] =
					lp->n_left->n_val.v_value.v_fval;
				p++;
			}
			lp = lp->n_right;
		}else if( lp->n_sym == INT_LIT ){
			tform[ p ] = lp->n_val.v_value.v_ival;
			p++;
			break;
		}else if( lp->n_sym == FLOAT_LIT ){
			tform[ p ] = lp->n_val.v_value.v_fval;
			p++;
			break;
		}
	}
	if( p != 6 )
		errormsg( TRUE, E_XFORM_PARM );
}
