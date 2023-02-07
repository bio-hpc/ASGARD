#define	G_UNDEF	(-1)

#define	G_REFERENCE	1

#define	G_PARALLEL	0
#define	G_ANTIPARALLEL	1

#define	G_NONE	0
#define	G_WATSONCRICK	1
#define	G_COPY	2

#define	G_START	0
#define	G_STOP	1	

typedef	struct	res_t	{
	struct	res_t	*r_next;
	char	*r_in;
	char	*r_out;
	char	*r_reslib;
	char	*r_resname;
	char	*r_reschar;
} RES_T;

typedef	struct	str_t	{
	char	*s_name;
	int	s_ref;
	TRANSFORM_T	*s_offset;
	int	s_orient;
	int	s_defres;
	int	s_n_res;
	RES_T	*s_res;
} STR_T;

typedef	struct	geom_t	{
	char	*g_name;
	TRANSFORM_T	g_helix;
	int	g_n_str;
	STR_T	*g_strs;
} GEOM_T;
