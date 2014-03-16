transporte-2d
=============
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "mpi.h"
//#include <Python.h>
//#include <iomanip>
using namespace std;
//funciones

double cond_inicial(double x0, double y0) {
	// Aqui podemos definir una funcion de x, lo de abajo es un parentesis logico si es cierto es 1 y lo demas es 0
	double ci;
	//ci=1.0*(y0<-0.5)+0.0*(y0>=-0.5);
	ci=1.0*(sqrt(x0*x0+y0*y0)<0.25);
	return ci;
}

double cond_contorno(double x0,double y0,double t) {
	// Aqui podemos definir una funcion para la cond. de contorno ,es la concentracion
	double c0;
	c0=0.0;
	return c0;
}

int main(int argc, char **argv) {/*es obligatorio que se llame main el programa principal, todolo que empieza por * es un puntero, ** son dos punteros, es un modelo dinámico*/
	double xmin,xmax,ymin,ymax,T,c1,c2,cfl,tiempo,x,y,t;/* es obligatorio el tipo*/
	int npx,npy,i,j,pij,pipj,pimj,pijm,pijp,nplx,nply,npl;
	char* fichero_salida;	/*cadena de caracteres*/
	double dx,dy,dt;
	FILE *fp;/*clase abstracta*/
	if (argc != 12) {/*esto es solo para tomar los datos,coge 8 donde el primero es el nombre del programa*/
		cerr<<"Uso:"<<endl;
		cerr<<argv[0] <<"xmin xmax ymin ymax T npx npy cfl c1 c2 fichero salida "<<endl;
		cerr<<"xmin: Comienzo del intervalo x."<<endl;
		cerr<<"xmax: Final del inervalo x."<<endl;
		cerr<<"ymin: Comienzo del intervalo y."<<endl;
		cerr<<"ymax: Final del inervalo y."<<endl;
		cerr<<"T: Tiempo total de integraci?n."<<endl;
		cerr<<"npx: N? de particiones en x"<<endl;
		cerr<<"npy: N? de particiones en y"<<endl;
		cerr<<"cfl: Coef. estabilidad."<<endl;
		cerr<<"c1: Coef. vel. de transporte"<<endl;
		cerr<<"c2: Coef. vel. de transporte"<<endl;
		cerr<<"Nombre fichero de salida"<<endl;
		return -1;
	};
	xmin=atof(argv[1]);
	xmax=atof(argv[2]);
	ymin=atof(argv[3]);
	ymax=atof(argv[4]);
	T=atof(argv[5]);
	npx=atoi(argv[6]);
	npy=atoi(argv[7]);
	cfl=atof(argv[8]);
	c1=atof(argv[9]);
	c2=atof(argv[10]);
	fichero_salida=argv[11];
	//vector <double> sol0((npx+1)*(npy+1));
       // vector <double> sol1((npx+1)*(npy+1));
	double sol0[(npx+1)*(npy+1)];
	double sol1[(npx+1)*(npy+1)];
	if((myid==numprocs-1)&&(((npx+1)%numprocs)!=0)){
		nplx=((npx+1)/numprocs)+((npx+1)%numprocs);
	}
	else{
		nplx=(npx+1)/numprocs;
	}
	
	if((myid==numprocs-1)&&(((npy+1)%numprocs)!=0)){
		nply=((npy+1)/numprocs)+((npx+1)%numprocs);
	}
	else{
		nply=(npy+1)/numprocs;
	}
	
	startime=MPI_Wtime();
	
        dx=(xmax-xmin)/double(npx);
	dy=(ymax-ymin)/double(npy);
        dt=cfl*min(dx/fabs(c1)+1e-15,dy/fabs(c2)+1e-15);
	printf("c1: %15.8f\n",c1);
	printf("c2: %15.8f\n",c2);
	printf("Dt: %15.8f\n",dt);
	
	for (j=0; j<=(npy); j++) {
		for (i=0; i<=(npx); i++) {
			double x;
			double y;
			x=xmin+double(i)*dx;
			y=ymin+double(j)*dy;
			sol0[i+j*(npx+1)]=cond_inicial(x,y);
			}
		}
		
	for(i=0;i<npl+2;i++){
		sol1[i]=sol0[i];
		}
	tiempo=0.0; 	
	
	fp=fopen(fichero_salida,"wt");
	fprintf(fp,"%15.8f %15.8f %5i %15.8f %15.8f %5i",xmin,xmax,npx,ymin,ymax,npy);
	fprintf(fp,"\n");
	
	fprintf(fp,"%15.8f",tiempo);
	for (j=0; j<=(npy); j++) {
		for (i=0; i<=(npx); i++) {
		fprintf(fp,"%15.8f",sol0[i+j*(npx+1)]);//por filas, tiempo y solucion de ese tiempo
		}
	}
	fprintf(fp,"\n");
	
	
	while (tiempo<T) {
		//bucle en tiempo
		printf("Tiempo: %15.8f\n",tiempo+dt);
		for (j=1; j<(npy); j++) {
			for (i=1; i<(npx); i++) {
				//bucle en espacio
				pij=i+j*(npx+1);
				pipj=i+1+j*(npx+1);
				pimj=i-1+j*(npx+1);
				pijm=i+(j-1)*(npx+1);
				pijp=i+(j+1)*(npx+1);
				
				sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*sol0[pipj]-c1*sol0[pimj]-fabs(c1)*(sol0[pipj]-2.0*sol0[pij]+sol0[pimj]))
						-dt/(2.0*dy)*(c2*sol0[pijp]-c2*sol0[pijm]-fabs(c2)*(sol0[pijp]-2.0*sol0[pij]+sol0[pijm]));
			}
		}
			
			
		for(i=1;i<npx;i++){
				j=0;
				pij=i+j*(npx+1);
				pipj=i+1+j*(npx+1);
				pimj=i-1+j*(npx+1);
				pijm=i+(j-1)*(npx+1);
				pijp=i+(j+1)*(npx+1);
				
				
				x=xmin+i*dx;
				y=ymin;//+j que es 0
				sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*sol0[pipj]-c1*sol0[pimj]-fabs(c1)*(sol0[pipj]-2.0*sol0[pij]+sol0[pimj]))
						-dt/(2.0*dy)*(c2*sol0[pijp]-c2*cond_contorno(x,y,t)-fabs(c2)*(sol0[pijp]-2.0*sol0[pij]+cond_contorno(x,y,t)));
				
				j=npy;
				pij=i+j*(npx+1);
				pipj=i+1+j*(npx+1);
				pimj=i-1+j*(npx+1);
				pijm=i+(j-1)*(npx+1);
				pijp=i+(j+1)*(npx+1);
			
				
				y=ymax;
				sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*sol0[pipj]-c1*sol0[pimj]-fabs(c1)*(sol0[pipj]-2.0*sol0[pij]+sol0[pimj]))
						-dt/(2.0*dy)*(c2*cond_contorno(x,y,t)-c2*sol0[pijm]-fabs(c2)*(cond_contorno(x,y,t)-2.0*sol0[pij]+sol0[pijm]));
				
				}
		for(j=1;j<npy;j++){
			
				i=0;
				pij=i+j*(npx+1);
				pipj=i+1+j*(npx+1);
				pimj=i-1+j*(npx+1);
				pijm=i+(j-1)*(npx+1);
				pijp=i+(j+1)*(npx+1);
				
				
				x=xmin;//+i que es 0
				y=ymin+j*dy;
				sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*sol0[pipj]-c1*cond_contorno(x,y,t)-fabs(c1)*(sol0[pipj]-2.0*sol0[pij]+cond_contorno(x,y,t)))
						-dt/(2.0*dy)*(c2*sol0[pijp]-c2*sol0[pijm]-fabs(c2)*(sol0[pijp]-2.0*sol0[pij]+sol0[pijm]));
				
				
				i=npx;
				pij=i+j*(npx+1);
				pipj=i+1+j*(npx+1);
				pimj=i-1+j*(npx+1);
				pijm=i+(j-1)*(npx+1);
				pijp=i+(j+1)*(npx+1);
				x=xmax;
				
				sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*cond_contorno(x,y,t)-c1*sol0[pimj]-fabs(c1)*(cond_contorno(x,y,t)-2.0*sol0[pij]+sol0[pimj]))
						-dt/(2.0*dy)*(c2*sol0[pijp]-c2*sol0[pijm]-fabs(c2)*(sol0[pijp]-2.0*sol0[pij]+sol0[pijm]));
				
				}
			
		//Esquinas
		//sol1[0]
		i=0;
		j=0;
		x=xmin;
		y=ymin;
		sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*sol0[pipj]-c1*cond_contorno(x,y,t)-fabs(c1)*(sol0[pipj]-2.0*sol0[pij]+cond_contorno(x,y,t)))
						-dt/(2.0*dy)*(c2*sol0[pijp]-c2*cond_contorno(x,y,t)-fabs(c2)*(sol0[pijp]-2.0*sol0[pij]+cond_contorno(x,y,t)));
		//sol1[npx]
		i=npx;
		j=0;
		x=xmax;
		y=ymin;
		sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*cond_contorno(x,y,t)-c1*sol0[pimj]-fabs(c1)*(cond_contorno(x,y,t)-2.0*sol0[pij]+sol0[pimj]))
						-dt/(2.0*dy)*(c2*sol0[pijp]-c2*cond_contorno(x,y,t)-fabs(c2)*(sol0[pijp]-2.0*sol0[pij]+cond_contorno(x,y,t)));
		//sol1[npx*npy]
		i=0;
		j=npy;
		x=xmin;
		y=ymax;
		sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*sol0[pipj]-c1*cond_contorno(x,y,t)-fabs(c1)*(sol0[pipj]-2.0*sol0[pij]+cond_contorno(x,y,t)))
						-dt/(2.0*dy)*(c2*cond_contorno(x,y,t)-c2*sol0[pijm]-fabs(c2)*(cond_contorno(x,y,t)-2.0*sol0[pij]+sol0[pijm]));
				
		//sol1[npx+npx*npy]
		i=npx;
		j=npy;
		x=xmax;
		y=ymax;
		sol1[pij]=sol0[pij]-dt/(2.0*dx)*(c1*cond_contorno(x,y,t)-c1*sol0[pimj]-fabs(c1)*(cond_contorno(x,y,t)-2.0*sol0[pij]+sol0[pimj]))
						-dt/(2.0*dy)*(c2*cond_contorno(x,y,t)-c2*sol0[pijm]-fabs(c2)*(cond_contorno(x,y,t)-2.0*sol0[pij]+sol0[pijm]));
				
		
		for(i=0;i<npl+2;i++){
			sol0[i]=sol1[i];
			}
		tiempo=tiempo+dt;
		fprintf(fp,"%15.8f",tiempo);
		for (j=0; j<=(npy); j++) {
			for (i=0; i<=(npx); i++) {
			fprintf(fp,"%15.8f",sol0[i+j*(npx+1)]);//por filas, tiempo y solucion de ese tiempo
			}
		}
	fprintf(fp,"\n");
		}
	
	
	fclose(fp);//cierro el fichero
	
	
	
	
	
	MPI_Finalize();
}
