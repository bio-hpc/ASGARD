#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"
#include "y.tab.h"

extern	int	debug;

extern	FILE	*cg_cfp;

int	avsfunctype = T_UNDEF;
char	*avsmodname = NULL;
static	AVSINFO_T	*avsinfo = NULL;
NODE_T	*avslist = NULL;	/* use to generate AVS mod	*/

static	void	CG_genavswidgets( NODE_T * );
static	void	CG_gen_1_avswidget( NODE_T *pp );
static	int	get_1_avswidgetlims( NODE_T *, int,
	VALUE_T *, VALUE_T *, VALUE_T *);
static	void	CG_genavscompplist( NODE_T * );
static	void	CG_genavscomppdefs( NODE_T * );
static	void	CG_genavscomplvars( NODE_T * );
static	void	CG_genavspchanges( NODE_T * );
static	void	CG_gen_1_avspchange( NODE_T * );
static	void	CG_gencallnabfunc( NODE_T * );

void	saveAVSinfo( char info[] )
{
	char	*wp;
	char	*ap;
	char	tp[ 20 ];
	char	np[ 40 ];
	char	*dp;
	int	nb, wl;
	AVSINFO_T	*avsp;

	*tp = '\0';
	*np = '\0';
	dp = NULL;
	for( ap = &info[ 9 ]; *ap; ){
		nb = strspn( ap, " \t\n" );
		wp = &ap[ nb ];
		wl = strcspn( wp, " \t\n" );
		if( !*tp ){
			strncpy( tp, wp, wl );
			tp[ wl ] = '\0';
		}else if( !*np ){
			strncpy( np, wp, wl );
			np[ wl ] = '\0';
		}else{
			dp = wp;
			break;
		}
		ap += nb + wl;
	}

	avsp = ( AVSINFO_T * )malloc( sizeof( AVSINFO_T ) );
	if( avsp == NULL ){
		fprintf( stderr, "saveAVSinfo: can't alloc AVSINFO_T\n" );
		return;
	}
	avsp->a_next = NULL;
	if( !strcmp( tp, "parm" ) )
		avsp->a_type = AI_PARM;
	else if( !strcmp( tp, "func" ) )
		avsp->a_type = AI_FUNC;
	else
		avsp->a_type = AI_UNDEF;
	ap = ( char * )malloc( (strlen( np ) + 1 ) * sizeof(char) );
	if( ap == NULL ){
		fprintf( stderr, "saveAVSinfo: can't alloc avsp->a_name\n" );
		return;
	}
	strcpy( ap, np );
	avsp->a_name = ap;
	ap = ( char * )malloc( (strlen( dp ) + 1) * sizeof(char) );
	if( ap == NULL ){
		fprintf( stderr, "saveAVSinfo: can't alloc avsp->a_data\n" );
		return;
	}
	strcpy( ap, dp );
	avsp->a_data = ap;

	avsp->a_next = avsinfo;
	avsinfo = avsp;
}

void	CG_genavswrapper( void )
{

	fprintf( cg_cfp, "\n" );
	fprintf( cg_cfp, "AVSinit_modules()\n" );
	fprintf( cg_cfp, "{\n" );
	fprintf( cg_cfp, "int %s_desc();\n", avsmodname );
	fprintf( cg_cfp, "\n" );
	fprintf( cg_cfp, "AVSmodule_from_desc(%s_desc);\n", avsmodname );
	fprintf( cg_cfp, "}\n" );

	fprintf( cg_cfp, "%s_desc()\n", avsmodname );
	fprintf( cg_cfp, "{\n" );
	fprintf( cg_cfp, "int %s_compute();\n", avsmodname );
	fprintf( cg_cfp, "int p;\n" );
	fprintf( cg_cfp, "\n" );

	fprintf( cg_cfp, "AVSset_module_name(\"%s\",MODULE_DATA);\n",
		avsmodname );
	fprintf( cg_cfp, "p = AVScreate_output_port(\"AtomNames\",\n" );
	fprintf( cg_cfp, "\"field 1D 1-space 1-vector uniform byte\" );\n" );
	fprintf( cg_cfp, "p = AVScreate_output_port(\"AtomBonds\",\n" );
	fprintf( cg_cfp, "\"field 1D 1-space 2-vector uniform integer\" );\n" );
	fprintf( cg_cfp, "p = AVScreate_output_port(\"AtomXYZ\",\n" );
	fprintf( cg_cfp, "\"field 1D 1-space 3-vector uniform float\" );\n" );

	CG_genavswidgets( avslist );

	fprintf( cg_cfp, "AVSset_compute_proc(%s_compute);\n", avsmodname );
	fprintf( cg_cfp, "return(1);\n" );
	fprintf( cg_cfp, "}\n" );

	fprintf( cg_cfp, "%s_compute(AtomNames,AtomBonds,AtomXYZ", avsmodname );
	CG_genavscompplist( avslist );
	fprintf( cg_cfp, ")\n" );

	fprintf( cg_cfp, "AVSfield_char **AtomNames;\n" );
	fprintf( cg_cfp, "AVSfield_int **AtomBonds;\n" );
	fprintf( cg_cfp, "AVSfield_float **AtomXYZ;\n" );
	CG_genavscomppdefs( avslist );

	fprintf( cg_cfp, "{\n" );
	fprintf( cg_cfp, "int dims0[1];\n" );
	CG_genavscomplvars( avslist );

	fprintf( cg_cfp, "int pc_%s = 0;\n", avsmodname );
	fprintf( cg_cfp, "MOLECULE_T *m_%s;\n", avsmodname );
	fprintf( cg_cfp, "\n" );

	CG_genavspchanges( avslist );

	fprintf( cg_cfp, "if(!pc_%s)\n", avsmodname );
	fprintf( cg_cfp, "return(0);\n" );

	CG_gencallnabfunc( avslist );

	fprintf( cg_cfp, "if(!nab2avs(m_%s,AtomNames,AtomBonds,AtomXYZ))\n",
		avsmodname );
	fprintf( cg_cfp, "return(0);\n" );

	fprintf( cg_cfp, "return(1);\n" );
	fprintf( cg_cfp, "}\n" );
}

static	void	CG_genavswidgets( NODE_T *avslist )
{
	NODE_T	*plp, *pp;

	for( plp = avslist; plp; plp = plp->n_right ){
		pp = plp->n_left;
		CG_gen_1_avswidget( pp );
	}
}

static	void	CG_gen_1_avswidget( NODE_T *pp )
{
	NODE_T	*tp;
	NODE_T	*ip;
	int	type;
	VALUE_T	def, lo, hi;

	if( pp->n_sym == SYM_LIST )
		pp = pp->n_left;
	tp = pp->n_left;
	type = tp->n_val.v_value.v_ival;
	ip = pp->n_right;
	if( ip->n_sym == SYM_LIST )
		ip = ip->n_left;
	if( ip->n_sym == SYM_INDIRECT )
		ip = ip->n_right;
	if( ip->n_sym == SYM_INDIRECT )
		ip = ip->n_right;
	if( get_1_avswidgetlims( ip, type, &def, &lo, &hi ) )
		return;
	switch( type ){
	case T_INT :
		fprintf( cg_cfp,
			"p = AVSadd_parameter(\"%s\",\"integer\",%d,%d,%d);\n",
			ip->n_val.v_value.v_cval,
			def.v_value.v_ival,
			lo.v_value.v_ival,
			hi.v_value.v_ival );
/*
		fprintf( cg_cfp,
		"AVSadd_parameter_prop(p,\"immediate\",\"boolean\",1);\n" );
*/
		fprintf( cg_cfp, "AVSconnect_widget(p,\"idial\");\n" );
		break;
	case T_FLOAT :
		fprintf( cg_cfp,
		"p = AVSadd_float_parameter(\"%s\",%8.1f,%8.1f,%8.1f);\n",
			ip->n_val.v_value.v_cval,
			def.v_value.v_fval,
			lo.v_value.v_fval,
			hi.v_value.v_fval );
/*
		fprintf( cg_cfp,
		"AVSadd_parameter_prop(p,\"immediate\",\"boolean\",1);\n" );
*/
		fprintf( cg_cfp, "AVSconnect_widget(p,\"dial\");\n" );
		break;
	case T_STRING :
		fprintf( cg_cfp,
	"p = AVSadd_parameter(\"%s\",\"string\",\"%s\",NULL,NULL);\n",
			ip->n_val.v_value.v_cval,
			def.v_value.v_cval );
/*       used to have extra arguments:
			lo.v_value.v_cval,
			hi.v_value.v_cval );   */
		fprintf( cg_cfp,
			"AVSadd_parameter_prop(p,\"width\",\"integer\",8);\n" );
		fprintf( cg_cfp, "AVSconnect_widget(p,\"typein\");\n" );
		break;
	}
}

static	int	get_1_avswidgetlims( NODE_T *ip, int type,
	VALUE_T *def, VALUE_T *lo, VALUE_T *hi )
{
	int		i_def, i_lo, i_hi;
	float		f_def, f_lo, f_hi;
	char		c_def[ 256 ];
	char		*id;
	AVSINFO_T	*avsp;

	id = ip->n_val.v_value.v_cval;
	for( avsp = avsinfo; avsp; avsp = avsp->a_next ){
		if( avsp->a_type != AI_PARM )
			continue;
		else if( strcmp( avsp->a_name, id ) )
			continue;
		else
			break;
	}

	switch( type ){
	case T_INT :
		if( !avsp ){
			i_def = 0;
			i_lo = -10000;
			i_hi = 10000;
		}else if( sscanf( avsp->a_data, "%d %d %d",
			&i_def, &i_lo, &i_hi ) != 3 ){
			i_def = 0;
			i_lo = -10000;
			i_hi = 10000;
		}
		def->v_value.v_ival = i_def;
		lo->v_value.v_ival = i_lo;
		hi->v_value.v_ival = i_hi;
		return( 0 );
	case T_FLOAT :
		if( !avsp ){
			f_def = 0.0;
			f_lo = -10000.0;
			f_hi = 10000.0;
		}else if( sscanf( avsp->a_data, "%f %f %f",
			&f_def, &f_lo, &f_hi ) != 3 ){
			i_def = 0;
			i_lo = -10000;
			i_hi = 10000;
		}
		def->v_value.v_fval = f_def;
		lo->v_value.v_fval = f_lo;
		hi->v_value.v_fval = f_hi;
		return( 0 );
	case T_STRING :
		if( !avsp )
			def->v_value.v_cval = "";
		else if( sscanf( avsp->a_data, "%s", c_def ) != 1 )
			def->v_value.v_cval = "";
		else
			def->v_value.v_cval = c_def;
		lo->v_value.v_cval = "";
		hi->v_value.v_cval = "";
		return( 0 );
	default :
		fprintf( stderr,
		"only ints, floats and strings can be parms to avs module\n"
			);
		return( 1 );
	}
}

static	void	CG_genavscompplist( NODE_T *avslist )
{
	NODE_T	*plp, *pp;
	NODE_T	*ip;

	for( plp = avslist; plp; plp = plp->n_right ){
		pp = plp->n_left;
		if( pp->n_sym == SYM_DECL ){
			pp = pp->n_right;
			ip = pp->n_left;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			fprintf( cg_cfp, ",%s", ip->n_val.v_value.v_cval );
		}
	}
}

static	void	CG_genavscomppdefs( NODE_T *avslist )
{
	NODE_T	*plp, *pp;
	NODE_T	*tp;
	NODE_T	*ip;
	int	type;

	for( plp = avslist; plp; plp = plp->n_right ){
		pp = plp->n_left;
		if( pp->n_sym == SYM_DECL ){
			tp = pp->n_left;
			type = tp->n_val.v_value.v_ival;
			pp = pp->n_right;
			ip = pp->n_left;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			switch( type ){
			case T_INT :
				fprintf( cg_cfp, "int" );
				break;
			case T_FLOAT :
				fprintf( cg_cfp, "float*" );
				break;
			case T_STRING :
				fprintf( cg_cfp, "char*" );
				break;
			}
			fprintf( cg_cfp, " %s;\n", ip->n_val.v_value.v_cval );
		}
	}
}

static	void	CG_genavscomplvars( NODE_T *avslist )
{
	NODE_T	*plp, *pp;
	NODE_T	*tp;
	NODE_T	*ip;
	VALUE_T	def, lo, hi;
	int	type;

	for( plp = avslist; plp; plp = plp->n_right ){
		pp = plp->n_left;
		if( pp->n_sym == SYM_DECL ){
			tp = pp->n_left;
			type = tp->n_val.v_value.v_ival;
			pp = pp->n_right;
			ip = pp->n_left;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			if( type == T_FLOAT ){
				fprintf( cg_cfp, "double l_%s = *%s;\n",
					ip->n_val.v_value.v_cval,
					ip->n_val.v_value.v_cval );
			}else if( type == T_STRING ){
				get_1_avswidgetlims(ip, type, &def, &lo, &hi);
				fprintf( cg_cfp,
					"static char l_%s[256] = \"%s\";\n",
					ip->n_val.v_value.v_cval,
					def.v_value.v_cval );
				fprintf( cg_cfp, "static char *lp_%s = l_%s;\n",
					ip->n_val.v_value.v_cval,
					ip->n_val.v_value.v_cval );
			}
		}
	}
}

static	void	CG_genavspchanges( NODE_T *avslist )
{
	NODE_T	*plp, *pp;

	for( plp = avslist; plp; plp = plp->n_right ){
		pp = plp->n_left;
		CG_gen_1_avspchange( pp );
	}
}

static	void	CG_gen_1_avspchange( NODE_T *pp )
{
	NODE_T	*tp;
	NODE_T	*ip;
	char	*pn;
	int	type;

	if( pp->n_sym == SYM_DECL ){
		tp = pp->n_left;
		type = tp->n_val.v_value.v_ival;
		pp = pp->n_right;
		ip = pp->n_left;
		if( ip->n_sym == SYM_INDIRECT )
			ip = ip->n_right;
		if( ip->n_sym == SYM_INDIRECT )
			ip = ip->n_right;
		pn = ip->n_val.v_value.v_cval;
		fprintf( cg_cfp, "if(AVSparameter_changed(\"%s\")){\n", pn );
		if( type == T_STRING ){
			fprintf( cg_cfp, "strcpy(l_%s,%s);\n", pn, pn );
			fprintf( cg_cfp,
		"AVSmodify_parameter(\"%s\",AVS_VALUE,\"\",NULL,NULL);\n",
				pn );
		}
		fprintf( cg_cfp, "pc_%s = 1;\n", avsmodname );
		fprintf( cg_cfp, "}\n" );
	}
}

static	void	CG_gencallnabfunc( NODE_T *avslist )
{
	NODE_T	*plp, *pp;
	NODE_T	*tp, *ip;

	fprintf( cg_cfp, "if(!(" );
	fprintf( cg_cfp, "m_%s = AVS_%s(", avsmodname, avsmodname );
	for( plp = avslist; plp; plp = plp->n_right ){
		pp = plp->n_left;
		if( pp->n_sym == SYM_DECL ){
			tp = pp->n_left;
			pp = pp->n_right;
			ip = pp->n_left;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			if( ip->n_sym == SYM_INDIRECT )
				ip = ip->n_right;
			if( tp->n_val.v_value.v_ival == T_INT )
				fprintf( cg_cfp, "&" );
			else if( tp->n_val.v_value.v_ival == T_FLOAT )
				fprintf( cg_cfp, "&l_" );
			else if( tp->n_val.v_value.v_ival == T_STRING )
				fprintf( cg_cfp, "&lp_" );
			fprintf( cg_cfp, "%s", ip->n_val.v_value.v_cval );
			if( plp->n_right )
				fprintf( cg_cfp, "," );
		}
	}
	fprintf( cg_cfp, ")))\n" );
	fprintf( cg_cfp, "return(0);\n" );
}
