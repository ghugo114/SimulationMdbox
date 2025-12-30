#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include <termios.h>
#include <unistd.h>
#include <ncurses.h>
#include "nrutil.h"
#include "gnuplot_i.h" 
#define DT 0.1
#define G 1.0e3
#define ME 9.11
#define QE 1.6
#define FNAMESIZE 64
#define PI 3.14159265358979323846
#define BOXL 15.
#define RA 30.
#define RG 50.
#define RS 0.00465
//#define GM 1.32712440018e20*pow(UAM,3)*pow(YRSE,2) Es lo mismo
#define GM 39.4769264142519
/**---------------Opción  0 anterior otras unidades (No se usan)---------------------------------------------------------------------**/			
#define MP 1e-9
#define MT 0.001
#define RP 0.00153172431509222 //Asumiendo misma densidad que titan y masa MP

/**--------Masas de los planetas exteriores relativas al sol (mass in simulation units)--------------**/
#define MJ 9.54791938e-4 
#define MS 2.85885980e-4
#define MU 4.36624404e-5
#define MN 5.15138902e-5
/**--------------------------------------------------------------------------------------------------**/



#define C 3.e3

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double mass;
    double charge;
    double rad;
  char name[FNAMESIZE];
} Particle;

typedef struct{
  int n,nreg;
  int foldi;
  double b,c;
  double gam,fa,fb,w,ttot,k;
  int nboxes,nb,num,mol,g,gopt,opt;
} PAR;

typedef struct {
    double gxmin, gxmax;
    double gymin, gymax;
    double gzmin, gzmax;
} Gnup;


typedef struct {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
} Box;

Gnup g;
PAR par;
Particle *particles,*p_p;
int n;
double *masa,*rad;
int nreg; //recording frequency


/** Allocate and free square matrices **/

double **square_matrix(int n) {
  double **mat;
  int i;
  mat = malloc(n * sizeof(double *));
  for(i = 0; i < n; i++)
    mat[i] = malloc(n * sizeof(double));
  return mat;
}
void free_square_matrix(double **mat, int n) {
  int i;
  for(i = 0; i < n ; i++)
    free(mat[i]);
  free(mat);
}

double **nonsquare_matrix(int n,int m) {
  double **mat;
  int i;
  mat = malloc(m * sizeof(double *));
  for(i = 0; i < m; i++)
    mat[i] = malloc(n * sizeof(double));
  return mat;
}



/** Funcion que realiza el producto de el vector por la matriz**/

double *mat_vect(double *v, double **mat,int n) {
  double *vt;
  int i, j;
  vt = malloc(n*sizeof(double));
  if(i == j)
    for(i = 0; i < n; i++)
      vt[i] = 0;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      vt[i]+=mat[i][j]*v[j];
  return(vt);
}


double vect_dot_vect(double *v, double *w) {
  double vw;
  int i, j;
 vw=0.0;
    for(i = 1; i <= 3; i++)
      vw += v[i]*w[i];
   return(vw);
}

double *vect_cross_vect(double *v, double *w){
	double *vw;
	int i,j;
	vw=dvector(1,3);
	vw[1]=v[2]*w[3]-v[3]*w[2];
	vw[2]=-(v[1]*w[3]-v[3]*w[1]);
	vw[3]=v[1]*w[2]-v[2]*w[1];
	return(vw);
}


void initialize_particle(Particle *p, double x, double y, double z,
             double vx, double vy, double vz, const char *name, double mass, double charge, double rad) {
  p->x = x;
  p->y = y;
  p->z = z;
  p->vx = vx;
  p->vy = vy;
  p->vz = vz;
  p->mass = mass;
  p->charge = charge;
  p->rad=rad;
  if(name) strncpy(p->name, name, FNAMESIZE-1);
  p->name[FNAMESIZE-1] = '\0';
}
void armo_rk4_v(double *v, Particle *p,int n){
  int i;
  int dim=3*n;
  for(i=1;i<=n;i++){
    v[(3*i)-2]=p[i-1].x;
    v[(3*i)-2+(dim)]=p[i-1].vx;
    v[(3*i)-1]=p[i-1].y;
    v[(3*i)-1+(dim)]=p[i-1].vy;
    v[(3*i)]=p[i-1].z;
    v[(3*i)+(dim)]=p[i-1].vz;

  }

 }

void coord_i(double *v, Particle *p, int part){
int i;
int dim=3*part;
for (i = 1;i<=part;i++ ){
    p[i-1].x=v[(3*i)-2];
    p[i-1].y=v[(3*i)-1];
    p[i-1].z=v[(3*i)];
    p[i-1].vx=v[(3*i)-2+(dim)];
    p[i-1].vy=v[(3*i)-1+(dim)];
    p[i-1].vz=v[(3*i)+(dim)];
  }
 }

void cart_to_polar(Particle *p, Particle *p_p,int part){
 int i;
 double r,phi,th,dr,dphi,dth;
 double x,y,z,vx,vy,vz;
 double **dpdx,*u_p,*u_c;
 u_p = malloc(3*sizeof(double));
 u_c = malloc(3*sizeof(double));

 dpdx=square_matrix(3);

 for (i = 0;i<part;i++ ){
	x=p[i].x;
	y=p[i].y;
	z=p[i].z;
	r=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	th=acos(z/r);
	phi=atan(y/x);
	p_p[i].x=r;
	p_p[i].y=th;
	p_p[i].z=phi;

	dpdx[0][0]=sin(p_p[i].y)*cos(p_p[i].z);
 	dpdx[0][1]=sin(p_p[i].y)*sin(p_p[i].z);
 	dpdx[0][2]=cos(p_p[i].y);
 	dpdx[1][0]=p_p[i].x*cos(p_p[i].y)*cos(p_p[i].z);
 	dpdx[1][1]=p_p[i].x*cos(p_p[i].y)*sin(p_p[i].z);
 	dpdx[1][2]=-p_p[i].x*sin(p_p[i].z);
 	dpdx[2][0]=-p_p[i].x*sin(p_p[i].y)*sin(p_p[i].z);
 	dpdx[2][1]=p_p[i].x*sin(p_p[i].y)*cos(p_p[i].z);
 	dpdx[2][2]=0.0;

     	
	u_c[0]=p[i].vx;
 	u_c[1]=p[i].vy;
 	u_c[2]=p[i].vz;
 	u_p=mat_vect(u_c, dpdx,3);
 	p_p[i].vx=u_p[0];
 	p_p[i].vy=u_p[1];
 	p_p[i].vz=u_p[2];
	}

free_square_matrix(dpdx, 3);
free(u_p);
free(u_c);
}

void polar_to_cart(Particle *p, Particle *p_p,int part){
 int i;	
 double **dpdx,*u_p,*u_c;
 u_p = malloc(3*sizeof(double));
 u_c = malloc(3*sizeof(double));

 dpdx=square_matrix(3);
 for (i = 0;i<part;i++ ){
    p[i].x=p_p[i].x*sin(p_p[i].y)*cos(p_p[i].z);
    p[i].y=p_p[i].x*sin(p_p[i].y)*sin(p_p[i].z);
    p[i].z=p_p[i].x*cos(p_p[i].y);
 
 dpdx[0][0]=sin(p_p[i].y)*cos(p_p[i].z);
 dpdx[1][0]=sin(p_p[i].y)*sin(p_p[i].z);
 dpdx[2][0]=cos(p_p[i].y);
 dpdx[0][1]=p_p[i].x*cos(p_p[i].y)*cos(p_p[i].z);
 dpdx[1][1]=p_p[i].x*cos(p_p[i].y)*sin(p_p[i].z);
 dpdx[2][1]=-p_p[i].x*sin(p_p[i].z);
 dpdx[0][2]=-p_p[i].x*sin(p_p[i].y)*sin(p_p[i].z);
 dpdx[1][2]=p_p[i].x*sin(p_p[i].y)*cos(p_p[i].z);
 dpdx[2][2]=0.0;
 u_p[0]=p_p[i].vx;
 u_p[1]=p_p[i].vy;
 u_p[2]=p_p[i].vz;
 u_c=mat_vect(u_p, dpdx,3);
 p[i].vx=u_c[0];
 p[i].vy=u_c[1];
 p[i].vz=u_c[2];
}
free_square_matrix(dpdx, 3);
free(u_p);
free(u_c);
}

void elem_orbitales(Particle *sis){
  FILE *arch;
  int i,exp,j,n,nplan;
  double t;
  
  double aa[8];

  if((arch=fopen("datos.dat", "r")) == NULL){
    printf("No se puede leer el archivo/%s\n","sist_solar.dat");
  }
  else
    printf(" - Datos de elementos orbitales cargados\n");
  while(fscanf(arch,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&i,&aa[0],&aa[1],&aa[2],&aa[3],&aa[4],&aa[5],&aa[6],&aa[7])!=EOF){
        
   
    printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,aa[0],aa[1],aa[2],aa[3],aa[4],aa[5],aa[6],aa[7]);
    sis[i].x = aa[0];
    sis[i].y = aa[1];
    sis[i].z = aa[2];       
    sis[i].vx = aa[3];
    sis[i].vy = aa[4];
    sis[i].vz = aa[5];
    sis[i].mass = aa[6];
    sis[i].rad = aa[7];   
 }
  fclose(arch);

}

void derivs(double x, double v[], double dv[]){
 int i,j,dim;  
 dim=par.n*3;
 double omega;
 omega=par.w;
 double *dvout;
 dvout=dvector(1,2*dim);
 int n =par.n;


for(i=1;i<=n;i++){   
	 if((v[(3*i)-2] <= RS)){ 
	  dvout[(3*i)-2+dim]=0.0;
	  dvout[(3*i)-2]=0.0;
  	  v[(3*i)-2]=RS;	 
	  dvout[(3*i)-1+dim]=0.0;
	
  	  dvout[(3*i)+dim]=0.0;
  	i++;
	 }
	 else if (v[(3*i)-2]>1000.){
	  dvout[(3*i)-2+dim]=0.0;
	  dvout[(3*i)-2]=0.0;
  	  v[(3*i)-2]=RS;	 
	  dvout[(3*i)-1+dim]=0.0;
	
  	  dvout[(3*i)+dim]=0.0;
  	i++;
	 }

	 
	 dvout[(3*i)-2+dim]=v[(3*i)-2]*(pow(v[(3*i)-1+dim],2)+(sin(v[(3*i)-1])*pow(v[(3*i)+dim],2))+(pow(sin(v[(3*i)-1]),2)*omega*(omega+v[(i*3)+dim])))-(GM/pow(v[(3*i)-2],2))-(par.b*par.gam*v[(3*i)-2+dim]);    
	 dvout[(3*i)-2]=v[(3*i)-2+(dim)];
	 dv[(3*i)-2]=dvout[(3*i)-2];
	 dv[(3*i)-2+dim]=dvout[(3*i)-2+dim];
	 
	 dvout[(3*i)-1+(dim)]=(sin(v[(3*i)-1])*cos(v[(3*i)-1])*(pow(v[(3*i)+dim]+omega,2)))-((2.*v[(3*i)-2+dim]*v[(3*i)-1+dim])/v[(3*i)-2])-(par.c*par.gam*v[(3*i)-1+dim]*v[(3*i)-2]);
	 dvout[(3*i)-1]=v[(3*i)-1+(dim)];
	 dv[(3*i)-1]=dvout[(3*i)-1];
	 dv[(3*i)-1+dim]=dvout[(3*i)-1+dim];
    
	 dvout[(3*i)+(dim)]=-2.*(omega+v[(3*i)+dim])*((v[(3*i)-1+dim]/tan(v[(3*i)-1]))+(v[(3*i)-2+dim]/v[(3*i)-2]))-(par.gam*par.c*v[(3*i)+dim]*v[(3*i)-2]*sin(v[(3*i)-1]));
   
	 dvout[(3*i)]=v[(3*i)+(dim)];
	 dv[(3*i)]=dvout[(3*i)];
	 dv[(3*i)+dim]=dvout[(3*i)+dim];
}


 free_dvector(dvout,1,dim*2);  
}






void rkdumb(double vstart[], int nvar, int nstep,
  void (*derivs)(double, double [], double [])){
  void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
  void (*derivs)(double, double [], double []));
  int i, j,count,nreg;
  long k;

  double d,x,h,fase,r12,r32;
  double *v,*vout,*dv;
  int n=par.n;
  int dim=par.n*3;
  char fname[FNAMESIZE];
  FILE *arch;
  arch = malloc(n * sizeof(FILE *));

  double *cenrgy,ue,ce;
 
  h=DT;
  x = 0.0; /* initialize time */
 count=1;
 

  v=dvector(1,nvar);
  vout=dvector(1,nvar);
  dv=dvector(1,nvar);
  for (i=1;i<=nvar;i++) {
    v[i]=vstart[i];
  }
  
  FILE *gnuplot = popen("gnuplot", "w");
  int paused = 0;
  struct termios orig_term;
  int orig_flags = -1;
 
  /* Setup gnuplot mouse support (optional) */
  fprintf(gnuplot, "set mouse\n");
  fflush(gnuplot);
  int l=0;
//  sprintf(fname,"particle_%d.dat",l);
 // arch=fopen(fname,"w");
 
  /* Make stdin non-canonical and non-blocking so we can toggle pause with keys */
  if (tcgetattr(STDIN_FILENO, &orig_term) == 0) {
    struct termios raw = orig_term;
    raw.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &raw);
    orig_flags = fcntl(STDIN_FILENO, F_GETFL, 0);
    fcntl(STDIN_FILENO, F_SETFL, orig_flags | O_NONBLOCK);

  }  

  h=DT;
  x = 0.0; /* initialize time */
 count=1;


  j=1;
 	g.gxmax=RG;
	g.gzmax=RG;
double r2,rij;
   for (k=1;k<nstep/100;k++){ 	

 	(*derivs)(x,v,dv);
    	rk4(v,dv,nvar,x,h,vout,derivs);
  
     if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
    x += h;
 
    for (i=1;i<=nvar;i++) {
      v[i]=vout[i];
    }

    coord_i(v,p_p,par.n);
    polar_to_cart(particles,p_p,par.n);
    
 g.gxmax=RG; 
 g.gzmax=RG;    
    fprintf(gnuplot, "set title 'Particles (x,y,z) at t=%lf' \n", x);
    fprintf(gnuplot, "set xrange [%lf:%lf]\nset yrange [%lf:%lf]\nset zrange [%lf:%lf]\n",-g.gxmax,g.gxmax,-g.gxmax,g.gxmax,-g.gzmax,g.gzmax);
    fprintf(gnuplot, "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\n");
    fprintf(gnuplot, "set isotropic\n");
    fprintf(gnuplot, "splot '-' with p pt 7 ps 0.5\n");
    for (i = 1; i <=n; i++) {
	    fprintf(gnuplot,"%lf\t%lf\t%lf\n", particles[i].x,particles[i].y,particles[i].z);	
    }
    fprintf(gnuplot, "e\n");
 
    fflush(gnuplot);
   
      /* Check keyboard for pause/resume/quit (non-blocking) */
    
      char ch;
      ssize_t r = read(STDIN_FILENO, &ch, 1);
      if (r == 1) {
        if (ch == ' ' || ch == 'p' || ch == 'P') {
          paused = !paused;
          if (paused) fprintf(stderr, "\n[PAUSED] Press 'p' or space to resume, 'q' to quit\n");
          else fprintf(stderr, "\nResuming...\n");
        } else if (ch == 'q' || ch == 'Q') {
              fprintf(stderr, "\nQuit requested. Exiting...\n");
              break;
            }
          }
         /* If paused, wait here until resumed or quit */
        while (paused) {
          /* Let gnuplot process events briefly */
          fprintf(gnuplot, "pause 0.01\n");
          fflush(gnuplot);
          char ch;
          ssize_t r = read(STDIN_FILENO, &ch, 1);
          if (r == 1) {
            if (ch == ' ' || ch == 'p' || ch == 'P') {
              paused = 0;
	      fprintf(stderr, "\nResuming...\n");
              break;
            } else if (ch == 'q' || ch == 'Q') {
              fprintf(stderr, "\nQuit requested. Exiting...\n");
             // goto cleanup_and_exit;
            }
          }
          usleep(100000);
        }
        
   }	
count=1;
for (k=(nstep/100)+1;k<95*nstep/100;k++){



	(*derivs)(x,v,dv);
    	rk4(v,dv,nvar,x,h,vout,derivs);
  
     if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
    x += h;
 
    for (i=1;i<=nvar;i++) {
      v[i]=vout[i];
    }
if(count==100){
	printf("%d\t%lf\n",k,x);
	count=0;
}
count++;
}

for (k=(95*nstep/100)+1;k<nstep;k++){ 	

	(*derivs)(x,v,dv);
    	rk4(v,dv,nvar,x,h,vout,derivs);
  
     if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
    x += h;
 
    for (i=1;i<=nvar;i++) {
      v[i]=vout[i];
    }

     coord_i(v,p_p,par.n);
    polar_to_cart(particles,p_p,par.n);

g.gxmax=RG/3; 
 g.gzmax=RG/3;;    
    fprintf(gnuplot, "set title 'Particles (x,y,z) at t=%lf' \n", x);
    fprintf(gnuplot, "set xrange [%lf:%lf]\nset yrange [%lf:%lf]\nset zrange [%lf:%lf]\n",-g.gxmax,g.gxmax,-g.gxmax,g.gxmax,-g.gzmax,g.gzmax);
    fprintf(gnuplot, "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\n");
   // fprintf(gnuplot, "set isotropic\n");
    fprintf(gnuplot, "splot '-' with p pt 7 ps 0.5\n");
    for (i = 0; i <n; i++) {
	    fprintf(gnuplot,"%lf\t%lf\t%lf\n", particles[i].x,particles[i].y,particles[i].z);	
    }
    fprintf(gnuplot, "e\n");
 
    fflush(gnuplot);
   
      /* Check keyboard for pause/resume/quit (non-blocking) */
    
      char ch;
      ssize_t r = read(STDIN_FILENO, &ch, 1);
      if (r == 1) {
        if (ch == ' ' || ch == 'p' || ch == 'P') {
          paused = !paused;
          if (paused) fprintf(stderr, "\n[PAUSED] Press 'p' or space to resume, 'q' to quit\n");
          else fprintf(stderr, "\nResuming...\n");
        } else if (ch == 'q' || ch == 'Q') {
              fprintf(stderr, "\nQuit requested. Exiting...\n");
              break;
            }
          }
         /* If paused, wait here until resumed or quit */
        while (paused) {
          /* Let gnuplot process events briefly */
          fprintf(gnuplot, "pause 0.01\n");
          fflush(gnuplot);
          char ch;
          ssize_t r = read(STDIN_FILENO, &ch, 1);
          if (r == 1) {
            if (ch == ' ' || ch == 'p' || ch == 'P') {
              paused = 0;
	      fprintf(stderr, "\nResuming...\n");
              break;
            } else if (ch == 'q' || ch == 'Q') {
              fprintf(stderr, "\nQuit requested. Exiting...\n");
              goto cleanup_and_exit;
            }
          }
          usleep(100000);
        }
        
   }

    coord_i(v,p_p,par.n);
    polar_to_cart(particles,p_p,par.n);


arch=fopen("datos.dat","w");
for(i=0;i<par.n;i++){
	fprintf(arch,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,particles[i].x,particles[i].y,particles[i].z,particles[i].vx,particles[i].vy,particles[i].vz,masa[i]*1.e9,rad[i]);
}
fclose(arch);







cleanup_and_exit:
  /* Restore terminal settings if we changed them */
  if (orig_flags != -1) {
    tcsetattr(STDIN_FILENO, TCSANOW, &orig_term);
    fcntl(STDIN_FILENO, F_SETFL, orig_flags);
  }
 
  if (gnuplot) pclose(gnuplot);
  if (arch) fclose(arch);



}

int main(int argc, char *argv[]) {

    
   
    if (argc < 7) {
		       
    printf("Uso: %s <N> <ttot> <w> <k> <b> <c> <gam>\n", argv[0]);
    exit(1);
	}
 
  // Lee los parametros requeridos desde la linea de comandos.
    par.n = atoi(argv[1]);
    par.ttot = atof(argv[2]);
    par.w = atof(argv[3]);
    par.k=atof(argv[4]);
    par.b=atof(argv[5]);
    par.c=atof(argv[6]);
    par.gam=atof(argv[7]);
    /*
   par.gopt=atoi(argv[8]);
    par.g=0;//atoi(argv[9]);
    par.mol=1;    
 */
    long steps;
    //FILE *arch1;   
    masa=dvector(1,par.n);
    rad=dvector(1,par.n);
    
    steps = par.ttot/DT;
    int rstep=0;
    par.nboxes=1;
    int n=par.n;
    int dim=6*par.n;
    particles = (Particle *)malloc(n * sizeof(Particle)); 
    p_p = (Particle *)malloc(n * sizeof(Particle)); 

    int opc=1;
    int i;
    double *vstart;
    vstart=dvector(1,dim);
    double th;     
        
    if(opc==0){
	elem_orbitales(p_p);
	cart_to_polar(p_p,particles,par.n);	
    }
    else{
	    for (i = 0; i < n/4; i++) {
	     th=(rand() / (double)RAND_MAX)*(int)(PI);
             
	     initialize_particle(&particles[i], 
                   RA,
		   th,
		   (rand() /  (double)RAND_MAX )*(int)(2.*PI),
                   0.0,
		   par.k*((rand() / (double)RAND_MAX - 0.5) * 2)/RA,
                   par.k*((rand() / (double)RAND_MAX - 0.5) * 2)/(RA*sin(th)),
                   "O",  MP,QE,RP);
	    }

    
	     for(i=1;i<=n;i++){
		     masa[i]=particles[i].mass;
		     rad[i]=particles[i].rad;
		        }
	     
    }	     
	     armo_rk4_v(vstart, particles, n);
             rkdumb(vstart,dim,steps,derivs);
 /* After simulation, show parameters again for editing before re-run */
    printf("\n\nSimulation completed. Edit parameters for next run?\n");
   // opc = 0;
      // free_dvector(masa,1,par.n);  	              	     
 //	free_dvector(rad,1,par.n);
  //   	free_dvector(vstart,1,dim);
      	return(0);
}
