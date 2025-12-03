#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include <termios.h>
#include <unistd.h>
#include "nrutil.h"
#include "gnuplot_i.h" 
#define DT 0.01
#define FNAMESIZE 64
#define PI 3.14159265358979323846

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double mass;
} Particle;
typedef struct{
  int n;
  int fold;
  double f,w,ttot,b,kc,kd,d,d_fold;
  double th,phi;
  int nboxes,nb,num;
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
Particle *particles;
int n;
double *masa;
double **conect;

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
void initialize_particle(Particle *p, double x, double y, double z, 
                         double vx, double vy, double vz, double mass) {
    p->x = x;
    p->y = y;
    p->z = z;
    p->vx = vx;
    p->vy = vy;
    p->vz = vz;
    p->mass = mass;
}
void armodexy(double *v, Particle *p,int part){
  int i;
  int dim=3*part;
  for(i=1;i<=part;i++){
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


void apply_boundary_conditions(Particle *p){
    for(int i=0;i<par.n;i++){
    if (p[i].x < 0.0 || p[i].x > 10.0) p[i].vx *= -1.;
    if (p[i].y < 0.0 || p[i].y > 10.0) p[i].vy *= -1.;
    if (p[i].z < 0.0 || p[i].z > 10.0) p[i].vz *= -1.;

    if (p[i].x < 0.00) p[i].x = 0.0;
    if (p[i].x > 10.0) p[i].x = 10.0;
    if (p[i].y < 0.00) p[i].y = 0.0;
    if (p[i].y > 10.0) p[i].y = 10.0;
    if (p[i].z < 0.00) p[i].z = 0.0;
    if (p[i].z > 10.0) p[i].z = 10.0;
    }
}

void record_positions(double t,Particle *particles, int n, FILE *arch) {
	for (int i = 0; i < n; i++){
	fprintf(arch, "%d\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%f\n",i,particles[i].x, particles[i].y, particles[i].z,particles[i].vx, particles[i].vy, particles[i].vz,particles[i].mass);
}
fclose(arch);
} 


    
void derivs(double x, double v[], double dv[]){
  int i,dim;  
  double roo,d,tau,gam,f,alfa;
  double *dvout;  
  double *fa,*fx,*fy,*fz;
  double *dd,*cc,*cx,*cy,*cz;
  dd=dvector(1,par.n+1);
  cc=dvector(1,par.n+1);
  int j;
  d=par.d;
  dim=3*par.n;  
  fa=dvector(1,3);
  fx=dvector(0,par.n);
  fy=dvector(0,par.n);
  fz=dvector(0,par.n);
  cx=dvector(1,par.n+1);
  cy=dvector(1,par.n+1);
  cz=dvector(1,par.n+1);
  fa[1]=sin(par.th)*cos(par.phi);
  fa[2]=sin(par.th)*sin(par.phi);
  fa[3]=cos(par.th);
  dvout=dvector(1,dim*2);  
  j = 1;

  /* Calculate forces to build linear conections of type:
        O=C=O=C=O=C=O...   */

      for(i=((j-1)*(par.n/par.nb))+1;i<=(par.n/par.nb)*j;i++){
      if((i==(par.n/par.nb)*j) && (j<par.nb)){
      j++;
      } 
      dd[i]=(1.-(d/(sqrt(pow((v[(3*i)-2]-v[(3*(i+1))-2]),2)+pow((v[(3*i)-1]-v[(3*(i+1))-1]),2)+pow((v[(3*i)]-v[(3*(i+1))]),2)))));
	    fx[i]=par.kd*dd[i]*(v[3*(i+1)-2]-v[(3*i)-2]);
      fy[i]=par.kd*dd[i]*(v[3*(i+1)-1]-v[(3*i)-1]);
      fz[i]=par.kd*dd[i]*(v[3*(i+1)]-v[(3*i)]);   
      
      /*Finish of linear conections*/
      /*Make bonds betwen segments ej:
      O=C=O=C=O=C
        | | | | |<-(this is par.fold=1 and this last bond is the one made before as =)          
      C=O=C=O=C=O
      */
      if(par.fold==1){
        if(i==((par.n/par.nb)*(j/2)) || i==((par.n/par.nb)*(j/2))-1) cc[i]=0.0; /*Fisrt and last are not conected*/
        else{
      cc[i]=(1.-((par.d_fold)/(sqrt(pow((v[(3*i)-2]-v[3*(j*(par.n/par.nb)-i+1)-2]),2)+pow((v[(3*i)-1]-v[((3*(j*(par.n/par.nb)-i+1))-1)]),2)+pow((v[(3*i)]-v[3*(j*(par.n/par.nb)-i+1)]),2)))));
      }
      cx[i]=par.kc*cc[i]*(v[(3*(j*(par.n/par.nb)-i+1))-2]-v[(3*i)-2]);
      cy[i]=par.kc*cc[i]*(v[(3*(j*(par.n/par.nb)-i+1))-1]-v[(3*i)-1]);
      cz[i]=par.kc*cc[i]*(v[(3*(j*(par.n/par.nb)-i+1))]-v[(3*i)]);
      }	
      cc[1]=0.0;
      cx[1]=0.0;
      cy[1]=0.0;
      cz[1]=0.0;
      cc[par.n]=0.0;
      cx[par.n]=0.0;
      cy[par.n]=0.0;
      cz[par.n]=0.0;

      if (par.fold==0){
      cc[i]=0.0;
      cx[i]=0.0;
      cy[i]=0.0;
      cz[i]=0.0;
      }
    }

      
      fx[0]=0.0;
      fy[0]=0.0;
      fz[0]=0.0;
      fz[(par.n)]=0.0;
      fx[(par.n)]=0.0;
      fy[(par.n)]=0.0;

    j=1; 
    

     for(i=(j-1)*(par.n/par.nb)+1;i<=(par.n/par.nb)*j;i++){
      if((i==((par.n/par.nb)*j))&&(j<par.nb)){
      j++;
      } 
    dvout[(3*i)-2+dim]=((fx[i]-fx[i-1]+cx[i])+fa[1]*par.f*cos(par.w*x)-((par.b*v[(3*i)-2+(dim)])))/masa[i];
    dvout[(3*i)-2]=v[(3*i)-2+(dim)];
    dv[(3*i)-2]=dvout[(3*i)-2];
    dv[(3*i)-2+dim]=dvout[(3*i)-2+dim];

    dvout[(3*i)-1+(dim)]=((fy[i]-fy[i-1]+cy[i])+fa[2]*par.f*cos(par.w*x)-((par.b*v[(3*i)-1+(dim)])))/masa[i] ;
    dvout[(3*i)-1]=v[(3*i)-1+(dim)];
    dv[(3*i)-1]=dvout[(3*i)-1];
    dv[(3*i)-1+dim]=dvout[(3*i)-1+dim];

    dvout[(3*i)+(dim)]=((fz[i]-fz[i-1]+cz[i])+fa[3]*par.f*cos(par.w*x)-((par.b*v[(3*i)+(dim)])))/masa[i];    dvout[(3*i)]=v[(3*i)+(dim)];
    dv[(3*i)]=dvout[(3*i)];
    dv[(3*i)+dim]=dvout[(3*i)+dim];
   
     }
      
  free_dvector(fa,1,3);
  free_dvector(fx,0,par.n+1);
  free_dvector(fy,0,par.n+1);
  free_dvector(fz,0,par.n+1);
    
  free_dvector(dd,1,par.n+1);
  free_dvector(cc,1,par.n+1);
  free_dvector(dvout,1,dim*2);  

  }

void rkdumb(double vstart[], int nvar, int nstep,
  void (*derivs)(double, double [], double [])){
  void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
  void (*derivs)(double, double [], double []));
  int i, j,count,registro;
  long k;

  
  double d,x,h,fase,r12,r32;
  double *v,*vx,*vy,*vout,*dv;
  double *cenrgy,ue,ce;
  char fname[FNAMESIZE];
 
     
  h=DT;
  x = 0.0; /* initialize time */
  FILE *arch;
  int n=par.n;
//  arch = malloc(n * sizeof(FILE *));
  FILE *gnuplot = popen("gnuplot", "w");
  int paused = 0;
  struct termios orig_term;
  int orig_flags = -1;

  /* Setup gnuplot mouse support (optional) */
  fprintf(gnuplot, "set mouse\n");
  fflush(gnuplot);

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

 // for (i=0;i<n;i++){  
//  sprintf(fname,"particle_%d.dat",i);
  arch=fopen("datos.dat","w");
 // }
  count=1;


  v=dvector(1,nvar);
  vout=dvector(1,nvar);
  dv=dvector(1,nvar);
  for (i=1;i<=nvar;i++) {
    v[i]=vstart[i];
  }


  for (k=1;k<nstep;k++){ 	
	
	  
	(*derivs)(x,v,dv);
    	rk4(v,dv,nvar,x,h,vout,derivs);
   
     if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
    x += h;
 
    for (i=1;i<=nvar;i++) {
      v[i]=vout[i];
      }
    coord_i(v, particles, n);
    apply_boundary_conditions(particles);
    armodexy(v,particles,n);
    
    /* Update axis ranges based on current particle positions */
    double xmin = particles[0].x, xmax = particles[0].x;
    double ymin = particles[0].y, ymax = particles[0].y;
    double zmin = particles[0].z, zmax = particles[0].z;
    for (i = 1; i < n; i++) {
      if (particles[i].x < xmin) xmin = particles[i].x;
      if (particles[i].x > xmax) xmax = particles[i].x;
      if (particles[i].y < ymin) ymin = particles[i].y;
      if (particles[i].y > ymax) ymax = particles[i].y;
      if (particles[i].z < zmin) zmin = particles[i].z;
      if (particles[i].z > zmax) zmax = particles[i].z;
    }
    /* Add padding to ranges */
    double xpad = (xmax - xmin) * 0.1;
    double ypad = (ymax - ymin) * 0.1;
    double zpad = (zmax - zmin) * 0.1;
    g.gxmin = xmin - xpad;
    g.gxmax = xmax + xpad;
    g.gymin = ymin - ypad;
    g.gymax = ymax + ypad;
    g.gzmin = zmin - zpad;
    g.gzmax = zmax + zpad;
    
 /* Send current x-y-z particle positions to gnuplot for live 3D plotting */
    fprintf(gnuplot, "set title 'Particles (x,y,z) at t=%lf %s (left-click to pause)'\n", x, paused ? "[PAUSED]" : "");
    fprintf(gnuplot, "set xrange [%lf:%lf]\nset yrange [%lf:%lf]\nset zrange [%lf:%lf]\n",g.gxmin,g.gxmax,g.gymin,g.gymax,g.gzmin,g.gzmax);
    fprintf(gnuplot, "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\n");
    //fprintf(gnuplot, "set view 60, 30\n");
   if(par.num==1) fprintf(gnuplot, "splot '-' with p pt 7 ps 1.5, '-' with labels notitle\n");
    else fprintf(gnuplot, "splot '-' with p pt 7 ps 1\n");
    for (i = 0; i < n; i++) {
      fprintf(gnuplot, "%lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].z);
    }
    fprintf(gnuplot, "e\n");

    if (par.num==1){
    /* Send particle indices as labels */
        for (i = 0; i < n; i++) {
      fprintf(gnuplot, "%lf %lf %lf %d\n", particles[i].x, particles[i].y, particles[i].z, i+1);
    }
    fprintf(gnuplot, "e\n");
    }
    fflush(gnuplot);

    /* Check keyboard for pause/resume/quit (non-blocking) */
    {
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
        count++; 

   
   
  }

  
  for (i = 0; i < n; i++){
    fprintf(arch,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,particles[i].x, particles[i].y, particles[i].z, particles[i].vx, particles[i].vy, particles[i].vz,particles[i].mass);
    
    }

cleanup_and_exit:
  /* Restore terminal settings if we changed them */
  if (orig_flags != -1) {
    tcsetattr(STDIN_FILENO, TCSANOW, &orig_term);
    fcntl(STDIN_FILENO, F_SETFL, orig_flags);
  }

  if (gnuplot) pclose(gnuplot);
  if (arch) fclose(arch);
  free_dvector(dv,1,nvar);
  free_dvector(vout,1,nvar);
  free_dvector(v,1,nvar);

}

  
/**/

int main(int argc, char *argv[]) {

    
   
    if (argc < 14) {
		       
    printf("Uso: %s <N> <ttot> <f_ext> <w_ext> <k(bond)> <k(fold)> <b(damp)> <d> <d(fold)> <th> <phi> <nb> <fold> <num_part>\n", argv[0]);
    exit(1);
	}

  // Lee los parametros requeridos desde la linea de comandos.
    par.n = atoi(argv[1]);
    par.ttot = atof(argv[2]);
    par.f=atof(argv[3]);
    par.w=atof(argv[4]);
    par.kd=atof(argv[5]);
    par.kc=atof(argv[6]);
    par.b=atof(argv[7]);
    par.d=atof(argv[8]);
    par.d_fold=atof(argv[9]);
    par.th=atof(argv[10]);
    par.phi=atof(argv[11]);
    par.nb=atoi(argv[12]);
    par.fold=atoi(argv[13]);
    par.num=atoi(argv[14]);
    //par.kc=par.kd/2; //constant for folded connections
    int steps;
        
    masa=dvector(1,par.n);
    steps = par.ttot/DT;
    int reg=1; //recording frequency
    int rstep=0;
    par.nboxes=1;
    int n=par.n;
    int n_particles = par.n;
    int dim=6*n_particles;
    particles = (Particle *)malloc(n * sizeof(Particle)); 
    g.gxmin=0.0;
    g.gxmax=10.0;
    g.gymin=0.0;
    g.gymax=10.0;
    g.gzmin=0.0;
    g.gzmax=10.0;
    int i;      

   // Initialize particles with random positions and velocities
    
    for (i = 0; i < n_particles; i++) {
        initialize_particle(&particles[i], 
                           rand() % 10, rand() % 10, rand() % 10,
                           (rand() / (double)RAND_MAX - 0.5) * 2,
                           (rand() / (double)RAND_MAX - 0.5) * 2,
                           (rand() / (double)RAND_MAX - 0.5) * 2,
                           ((1+pow(-1,i))/2)*2.66+((1+pow(-1,i+1))/2)*1.99);
    }

  

    for (i = 1; i <=n; i++){
        masa[i]=particles[i-1].mass;
    } 
    double *vstart, *vcont;
    double nga,nf,nw,nth,nphi,gxmax,gxmin,gymax,gymin,gzmax,gzmin;
    vstart=dvector(1,dim);
    vcont=dvector(1,dim);
    
       
    armodexy(vstart, particles, n);
    int opc,gopc;
    rkdumb(vstart, dim, steps, derivs);
 /*
    opc=1;
  scanf("Sigo igual?",&opc);
  while(opc==1){
  printf("Continúo con las coordenadas actuales cambiando parametros o leo sistema.dat?(1 sigo - 2 leo cord): ");
  scanf("%d", &opc);
  printf("Parámetros anteriores - ga,f_ext,w_ext: %lf\t%lf\t%lf\t%lf\t%lf\n",par.b,par.f,par.w,par.th,par.phi);
  scanf("%lf,%lf,%lf,%lf,%lf",&nga,&nf,&nw,&nth,&nphi);
  printf("Nuevos parámetros ingresados continúo %lf\n: ",par.ttot);
  par.b=nga;
  par.f=nf;
  par.w=nw;
  par.th=nth;
  par.phi=nphi;
  
  armodexy(vcont, particles, n);
  rkdumb(vcont, dim, steps, derivs);

  
  }*/

 free_dvector(masa,1,par.n);
 return 0;
}
  
