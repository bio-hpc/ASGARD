int removetranrot(x,v,minv)
REAL_T *x,*v,*minv;

{
  int n,i,j,nth;
  REAL_T Lx,Ly,Lz;
  REAL_T xs,ys,zs;
  REAL_T vsx,vsy,vsz;
  REAL_T xx,xy,xz,yy,yz,zz;
  REAL_T I[9];
  REAL_T Omx,Omy,Omz;
  REAL_T det,cofac;
  static REAL_T m_t;
  #define SCFAC 1.e-3
  n=3*prm->Natom;
  /*reset COM*/
  xs=ys=zs=vsx=vsy=vsz=0.0;

  if( m_t == 0.0)
    {
      for(i=0;i<n/3;i++)
	{
	  m_t += 1./minv[3*i];
	}
    }


  for(i=0;i<n/3;i++)
    {
      xs += x[3*i  ] / minv[3*i  ];
      ys += x[3*i+1] / minv[3*i+1];
      zs += x[3*i+2] / minv[3*i+2];
      vsx += v[3*i  ] / minv[3*i  ];
      vsy += v[3*i+1] / minv[3*i+1];
      vsz += v[3*i+2] / minv[3*i+2];
      }
  xs /= m_t;
  ys /= m_t;
  zs /= m_t;
  vsx /= m_t;
  vsy /= m_t;
  vsz /= m_t;
  
   for(i=0;i<n/3;i++)
    {
      x[3*i  ] -= xs;
      x[3*i+1] -= ys;
      x[3*i+2] -= zs;

      v[3*i  ] -= vsx;
      v[3*i+1] -= vsy;
      v[3*i+2] -= vsz;
      
      
    }

   /*now COM translation has been removed, he he he */
   /* now comes the harder part */
   /* we'll remove rotation  about COM as well */

   /*Angular momentum first*/
   Lx=Ly=Lz=0.0;
   for(i=0;i<n/3;i++)
     {
       Lx += SCFAC*(x[3*i+1]*v[3*i+2] - x[3*i+2]*v[3*i+1])/minv[3*i];
       Ly += SCFAC*(x[3*i+2]*v[3*i  ] - x[3*i  ]*v[3*i+2])/minv[3*i];
       Lz += SCFAC*(x[3*i  ]*v[3*i+1] - x[3*i+1]*v[3*i  ])/minv[3*i];
     }
   /*then tensor of inertia*/

   xx=xy=xz=yy=yz=zz=0.0;

   for(i=0;i<n/3;i++)
     {
       xx += x[3*i  ]*x[3*i  ]/minv[3*i];
       xy += x[3*i  ]*x[3*i+1]/minv[3*i];
       xz += x[3*i  ]*x[3*i+2]/minv[3*i];
       yy += x[3*i+1]*x[3*i+1]/minv[3*i];
       yz += x[3*i+1]*x[3*i+2]/minv[3*i];
       zz += x[3*i+2]*x[3*i+2]/minv[3*i];
     }
   I[0] = (yy+zz)*SCFAC;
   I[1] = -xy*SCFAC;
   I[2] = -xz*SCFAC;
   I[3] = -xy*SCFAC;
   I[4] = (xx+zz)*SCFAC;
   I[5] = -yz*SCFAC;
   I[6] = -xz*SCFAC;
   I[7] = -yz*SCFAC;
   I[8] = (xx+yy)*SCFAC;


   psinv(I,3);

   Omx=Omy=Omz=0.0;

   Omx = I[0]*Lx + I[1]*Ly + I[2]*Lz;
   Omy = I[3]*Lx + I[4]*Ly + I[5]*Lz;
   Omz = I[6]*Lx + I[7]*Ly + I[8]*Lz;

   /* the molecules 'angular velocity'*/

   /*   printf("ABSOMEGA %e\n", sqrt(Omx*Omx+Omy*Omy+Omz*Omz));*/
for(i=0;i<n/3;i++)
    {
      v[3*i  ] = v[3*i  ] - Omy*x[3*i+2] + Omz*x[3*i+1];
      v[3*i+1] = v[3*i+1] - Omz*x[3*i ] + Omx*x[3*i+2];
      v[3*i+2] = v[3*i+2] - Omx*x[3*i+1] + Omy*x[3*i  ];
      
      
    }


 return (0);
}

int psinv(double *v,int n)
     /*from the ccmath library*/
{ double z,*p,*q,*r,*s,*t; int j,k;
 for(j=0,p=v; j<n ;++j,p+=n+1){
   for(q=v+j*n; q<p ;++q) *p-= *q* *q;
   if(*p<=0.) return -1;
   *p=sqrt(*p);
   for(k=j+1,q=p+n; k<n ;++k,q+=n){
     for(r=v+j*n,s=v+k*n,z=0.; r<p ;) z+= *r++ * *s++;
     *q-=z; *q/= *p;
   }
 }
