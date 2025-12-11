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
#define DT 0.05
#define Q 1.6
#define FNAMESIZE 64
#define PI 3.14159265358979323846

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double mass;
    double charge;
  char name[FNAMESIZE];
} Particle;
typedef struct{
  int n;
  int fold;
  double f,w,ttot,b,kc,kd,kang,d,d_fold;
  double th,phi;
  int nboxes,nb,num,mol;
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
             double vx, double vy, double vz, const char *name, double mass, double charge) {
  p->x = x;
  p->y = y;
  p->z = z;
  p->vx = vx;
  p->vy = vy;
  p->vz = vz;
  p->mass = mass;
  p->charge= charge;
  if(name) strncpy(p->name, name, FNAMESIZE-1);
  p->name[FNAMESIZE-1] = '\0';
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

/* Return simulation-unit mass for an element symbol.
   Adjust the values below to match your simulation units. */
double sim_mass_for(const char *symbol) {
  if (!symbol) return 1.0;
  if (strcmp(symbol, "C") == 0) return 1.99;    /* Carbon (sim units) */
  if (strcmp(symbol, "N") == 0) return 2.32;    /* Nitrogen (sim units) - example */
  if (strcmp(symbol, "O") == 0) return 2.66;    /* Oxygen (sim units) - example */
  if (strcmp(symbol, "H") == 0) return 0.14;    /* Hydrogen (sim units) - example */
  if (strcmp(symbol, "S") == 0) return 32.06;   /* Sulfur (placeholder) */
  return 1.0; /* default */
}

/* Build a simple repeating unit: C-N-C-C-N-C pattern.
   The function fills the provided particles array in-place for up to
   `residues * 6` atoms. Positions are placed along +X for simplicity.
*/
void build_chain(Particle *particles, int residues) {
  if (!particles || residues <= 0) return;
  const char *pattern[3] = {"N", "C", "C"};
  int q[3];
  q[0]=1;
  q[1]=-1;
  q[2]=0;
  int atoms_per_res = 3;
  int total_atoms = residues * atoms_per_res;
  double x = 5, y = 5, z = 5; /* centered in box */
  double spacing = 1.5; /* approximate bond spacing */
  int idx = 0;
  int count=0;
  for (int r = 0; r < residues; r++) {
    for (int j = 0; j < atoms_per_res; j++) {
      const char *sym = pattern[j];
      double mass = sim_mass_for(sym);
         initialize_particle(&particles[idx],
                     x, y, z ,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   sym, mass,q[j]);
          x += spacing; 
          count++;
          idx++;
    }
    x+= 0.4;
    z += 0.4;
   
  }
}

void build_co2(Particle *particles, int residues) {
  if (!particles || residues <= 0) return;
  const char *pattern[3] = {"O", "C", "O"};
  int q[3];
  q[0]=1;
  q[1]=-1;
  q[2]=0;
  int atoms_per_res = 3;
  int total_atoms = residues * atoms_per_res;
  double x = 5, y = 5, z = 5; /* centered in box */
  double spacing = 1.5; /* approximate bond spacing */
  int idx = 0;
  int count=0;
  for (int r = 0; r < residues; r++) {
    for (int j = 0; j < atoms_per_res; j++) {
      const char *sym = pattern[j];
      double mass = sim_mass_for(sym);
         initialize_particle(&particles[idx],
                     x, y, z ,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   sym, mass,q[j]);
          x += spacing; 
          count++;
          idx++;
    }
    x+= 0.4;
    z += 0.4;
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
  fprintf(arch, "%d\t%s\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%f\n", i, particles[i].name, particles[i].x, particles[i].y, particles[i].z, particles[i].vx, particles[i].vy, particles[i].vz, particles[i].mass);
}
fclose(arch);
} 

 
    
void derivs(double x, double v[], double dv[]){
  int i,dim;  
  double roo,d,tau,gam,f,alfa;
  double *dvout;  
  double *fa,*fx,*fy,*fz,*fangx,*fangy,*fangz;
  double *rr,*dd,*cc,*cx,*cy,*cz;
  dd=dvector(1,par.n+1);
  cc=dvector(1,par.n+1);
  rr=dvector(1,par.n+1);
  int j;
  d=par.d;
  double cos_t;
  dim=3*par.n;  
  fa=dvector(1,3);
  fx=dvector(0,par.n);
  fy=dvector(0,par.n);
  fz=dvector(0,par.n);
  fangx=dvector(0,par.n);
  fangy=dvector(0,par.n);
  fangz=dvector(0,par.n);

  cx=dvector(1,par.n+1);
  cy=dvector(1,par.n+1);
  cz=dvector(1,par.n+1);
  fa[1]=sin(par.th)*cos(par.phi);
  fa[2]=sin(par.th)*sin(par.phi);
  fa[3]=cos(par.th);
  dvout=dvector(1,dim*2);  
  j = 1;
  cos_t=0.866;
  /* Calculate forces to build linear conections of type:
        O=C=O=C=O=C=O...   */

      for(i=((j-1)*(par.n/par.nb))+1;i<=(par.n/par.nb)*j;i++){
      if((i==(par.n/par.nb)*j) && (j<par.nb)){
      j++;
      } 
      rr[i]=sqrt(pow((v[(3*i)-2]-v[(3*(i+1))-2]),2)+pow((v[(3*i)-1]-v[(3*(i+1))-1]),2)+pow((v[(3*i)]-v[(3*(i+1))]),2));
      dd[i]=(1.-(d/(rr[i])));
	    
      fx[i]=par.kd*dd[i]*(v[3*(i+1)-2]-v[(3*i)-2]);
      fy[i]=par.kd*dd[i]*(v[3*(i+1)-1]-v[(3*i)-1]);
      fz[i]=par.kd*dd[i]*(v[3*(i+1)]-v[(3*i)]);   
      fangx[i]=par.kang*(((v[(3*i+2)-2]-v[(3*(i+1))-2])/rr[i+1])-cos_t*(v[(3*i)-2]-v[(3*(i+1))-2])/(rr[i]));    
      fangy[i]=par.kang*(((v[(3*i+2)-1]-v[(3*(i+1))-1])/rr[i+1])-cos_t*(v[(3*i)-1]-v[(3*(i+1))-1])/(rr[i]));
      fangz[i]=par.kang*(((v[(3*i+2)]-v[(3*(i+1))])/rr[i+1])-cos_t*(v[(3*i)]-v[(3*(i+1))])/(rr[i]));
      /*Finish of linear conections*/
      /*Make bonds betwen segments ej:
      C=N=C=C=N=C
        | | | | |<-(this is par.fold=1 and this last bond is the one made before as =)          
                C=N=C=C=N=C
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
    dvout[(3*i)-2+dim]=((fx[i]-fx[i-1]+cx[i])+particles[i].charge*fa[1]*par.f*cos(par.w*x)-((par.b*v[(3*i)-2+(dim)])))/masa[i];
    dvout[(3*i)-2]=v[(3*i)-2+(dim)];
    dv[(3*i)-2]=dvout[(3*i)-2];
    dv[(3*i)-2+dim]=dvout[(3*i)-2+dim];

    dvout[(3*i)-1+(dim)]=((fy[i]-fy[i-1]+cy[i])+particles[i].charge*fa[2]*par.f*cos(par.w*x)-((par.b*v[(3*i)-1+(dim)])))/masa[i] ;
    dvout[(3*i)-1]=v[(3*i)-1+(dim)];
    dv[(3*i)-1]=dvout[(3*i)-1];
    dv[(3*i)-1+dim]=dvout[(3*i)-1+dim];

    dvout[(3*i)+(dim)]=((fz[i]-fz[i-1]+cz[i])+particles[i].charge*fa[3]*par.f*cos(par.w*x)-((par.b*v[(3*i)+(dim)])))/masa[i];    
    dvout[(3*i)]=v[(3*i)+(dim)];
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

/* Display and edit simulation parameters using ncurses */
void display_and_edit_parameters(PAR *par) {
  initscr();
  cbreak();
  noecho();
  keypad(stdscr, TRUE);

  int max_y, max_x;
  getmaxyx(stdscr, max_y, max_x);

  WINDOW *param_win = newwin(30, 60, 2, (max_x - 60) / 2);
  box(param_win, 0, 0);

  int selected = 0;
  int running = 1;

  while (running) {
    wclear(param_win);
    box(param_win, 0, 0);
    mvwprintw(param_win, 1, 2, "=== SIMULATION PARAMETERS ===");
    mvwprintw(param_win, 2, 2, "Use UP/DOWN arrows to select, ENTER to edit, 'q' to exit");
    
    int row = 4;
    char *field_names[] = {
      "n (particles)", "ttot (total time)", "f (ext force)", "w (frequency)",
      "k (bond const)", "k (folding","b (damping)", "d (chain)", "d_fold (folding)",
      "th (theta)", "phi (phi)", "nb (num boxes)", "fold (folding)","num(swow nums=1)"
    };
    
    mvwprintw(param_win, row++, 4, "n:        %d", par->n);
    mvwprintw(param_win, row++, 4, "ttot:     %.4f", par->ttot);
    mvwprintw(param_win, row++, 4, "f:        %.4f", par->f);
    mvwprintw(param_win, row++, 4, "w:        %.4f", par->w);
    mvwprintw(param_win, row++, 4, "kd:        %.4f", par->kd);
    mvwprintw(param_win, row++, 4, "kc:        %.4f", par->kc);
    mvwprintw(param_win, row++, 4, "b:        %.4f", par->b);
    mvwprintw(param_win, row++, 4, "d:      %.4f", par->d);
    mvwprintw(param_win, row++, 4, "d_fold:        %.4f", par->d_fold);
    mvwprintw(param_win, row++, 4, "th:       %.4f", par->th);
    mvwprintw(param_win, row++, 4, "phi:      %.4f", par->phi);
    mvwprintw(param_win, row++, 4, "nb:       %d", par->nb);
    mvwprintw(param_win, row++, 4, "fold:     %d", par->fold);
    mvwprintw(param_win, row++, 4, "num:      %d", par->num);

    /* Highlight selected line */
    int highlight_row = 4 + selected;
    mvwchgat(param_win, highlight_row, 4, 40, A_REVERSE, 0, NULL);

    wrefresh(param_win);

    int ch = getch();
    if (ch == KEY_UP && selected > 0) selected--;
    else if (ch == KEY_DOWN && selected < 13) selected++;
    else if (ch == '\n') {
      /* Edit selected parameter */
      char input[20];
      echo();
      mvwprintw(param_win, 28, 4, "Enter new value: ");
      wgetstr(param_win, input);
      noecho();

      if (strlen(input) > 0) {
        switch (selected) {
          case 0: par->n = atoi(input); break;
          case 1: par->ttot = atof(input); break;
          case 2: par->f = atof(input); break;
          case 3: par->w = atof(input); break;
          case 4: par->kd = atof(input); break;
          case 5: par->kc = atof(input); break;
          case 6: par->b = atof(input); break;
          case 7: par->d = atof(input); break;
          case 8: par->d_fold = atof(input); break;
          case 9: par->th = atof(input); break;
          case 10: par->phi = atof(input); break;
          case 11: par->nb = atoi(input); break;
          case 12: par->fold = atoi(input); break;
          case 13: par->num = atoi(input); break;
        }
      }
    } else if (ch == 'q' || ch == 'Q') {
      running = 0;
    }
  }

  delwin(param_win);
  endwin();
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
     /*  Update axis ranges based on current particle positions */
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
  /*     Add padding to ranges */
    double xpad = (xmax - xmin) * 0.1;
    double ypad = (ymax - ymin) * 0.1;
    double zpad = (zmax - zmin) * 0.1;
  
    g.gxmin = xmin - xpad;
    g.gxmax = xmax + xpad;
    g.gymin = ymin - ypad;
    g.gymax = ymax + ypad;
    g.gzmin = zmin - zpad;
    g.gzmax = zmax + zpad;
   }
   /* Send current x-y-z particle positions to gnuplot for live 3D plotting */
    fprintf(gnuplot, "set title 'Particles (x,y,z) at t=%lf %s (left-click to pause)'\n", x, paused ? "[PAUSED]" : "");
    fprintf(gnuplot, "set xrange [%lf:%lf]\nset yrange [%lf:%lf]\nset zrange [%lf:%lf]\n",g.gxmin,g.gxmax,g.gymin,g.gymax,g.gzmin,g.gzmax);
    fprintf(gnuplot, "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\n");
    //fprintf(gnuplot, "set view 60, 30\n");
    if(par.num==1) fprintf(gnuplot, "splot '-' with lp pt 7 ps 2, '-' with labels notitle\n");
    if(par.num==2) fprintf(gnuplot, "splot '-' with lp pt 7 ps 2, '-' with labels notitle\n");

    else fprintf(gnuplot, "splot '-' with lp pt 7 ps 2\n");
    for (i = 0; i < n; i++) {
      fprintf(gnuplot, "%lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].z);
    }
    fprintf(gnuplot, "e\n");

    if ((par.num==2) || (par.num==3)){
    /* Send particle indices as labels */
        for (i = 0; i < n; i++) {
      fprintf(gnuplot, "%lf %lf %lf %s\n", particles[i].x, particles[i].y, particles[i].z, particles[i].name);
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

    
   
    if (argc < 15) {
		       
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
    par.mol=atoi(argv[15]); //1 for polipeptide, 2 for CO2 chain, else random particles 
    //par.kc=par.kd/2; //constant for folded connections
    if(par.mol==1) printf(" - Building alanine-like polypeptide chain\n");
    else if(par.mol==2) {
      printf(" - Building CO2 chain\n");
      par.kang=par.kc; //angular constant
    }
    else printf(" - Building random particle system\n");  
    int steps;
    FILE *arch1;   
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
    g.gxmax=20.0;
    g.gymin=0.0;
    g.gymax=20.0;
    g.gzmin=0.0;
    g.gzmax=20.0;
    int i;      

     // Initialize particles: either build an alanine-like chain (if par.num==2)
     // or random positions/velocities otherwise.
    if (par.mol == 1) {
      int residues = par.n / 3; /* N-C-C pattern -> 3 atoms per residue */
      if (residues < 1) residues = 1;
      build_chain(particles, residues);
    } 
    else if (par.mol ==2) {
      int residues = par.n / 3; /* O-C-O pattern -> 3 atoms per residue */
      if (residues < 1) residues = 1;
      build_co2(particles, residues);
    } 
    else {
      for (i = 0; i < n_particles; i++) {
        double mass = ((1+pow(-1,i))/2)*1.99+((1+pow(-1,i+1))/2)*2.66;
        initialize_particle(&particles[i], 
                   rand() % 10, rand() % 10, rand() % 10,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   (rand() / (double)RAND_MAX - 0.5) * 2,
                   "C", mass,pow(-1,i));
      }
    }
 /*    Abrir Condiciones iniciales 
  if(par.num==-1){
    arch1=fopen("datos.dat","r");
    
    printf(" - Datos de condiciones corrida anterior y masas cargados\n");
    int k;
    double tfin,aa[6];  
  while(fscanf(arch1,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&k,&aa[1],&aa[2],&aa[3],&aa[4],&aa[5],&tfin)!=EOF){
      particles[k].x=aa[1];
      particles[k].y=aa[2];
      particles[k].z=aa[3];
      particles[k].vx=aa[4];
      particles[k].vy=aa[5];
      particles[k].vz=tfin;
      particles[k].mass=aa[6];
    }     
  }
  fclose(arch1);
  */
    for (i = 1; i <=n; i++){
        masa[i]=particles[i-1].mass;
    } 
    double *vstart, *vcont;
    double nga,nf,nw,nth,nphi,gxmax,gxmin,gymax,gymin,gzmax,gzmin;
    vstart=dvector(1,dim);
    vcont=dvector(1,dim);
    
       
    armodexy(vstart, particles, n);
    int opc,gopc;
    
    /* Display parameters on startup */
    display_and_edit_parameters(&par);
    
    rkdumb(vstart, dim, steps, derivs);
    
    /* After simulation, show parameters again for editing before re-run */
    printf("\n\nSimulation completed. Edit parameters for next run?\n");
    opc = 1;
    while(opc == 1) {
      display_and_edit_parameters(&par);
      printf("Run again? (1=yes, 0=no): ");
      scanf("%d", &opc);
      if (opc == 1) {
        steps = par.ttot / DT;
        armodexy(vstart, particles, n);
        rkdumb(vstart, dim, steps, derivs);
      }
    }

 free_dvector(masa,1,par.n);
 return 0;
}
  
