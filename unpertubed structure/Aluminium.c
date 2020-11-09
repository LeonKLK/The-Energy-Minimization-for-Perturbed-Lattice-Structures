//
//  main.c
//  3_part1
//
//  Created by Leon Kwok on 10/10/2019.
//  Copyright © 2019 Leon Kwok. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double vb[3][3], inter_cp[3];
double inter_dp, h;
double sigma, epsilon;
int main() {
    double a;   //lattice constant
    int i,j,k,NoA, N;
    int nx, ny, nz;
    double total_energy;
    void cross_product(double,double,double,double,double,double);
    void dot_product(double,double,double,double,double,double);
    double lennard_jones(double);
    
    //the calculation below is Aluminium
    sigma=2.6175;   //using angstroms as unit
    epsilon=4577.6*8.617333262145e-5;     //4577.6*boltzmann constant eV/K
    nx=5;
    ny=5;
    nz=5;
    a=4.046;    //using angstroms as unit

    NoA=nx*ny*nz*4;

    double p[NoA][3];
    N=0;
    for (k=0; k<nz; ++k) {
           for (j=0; j<ny; ++j) {
               for (i=0; i<nx; ++i) {
                   p[N][0]=i*a;
                   p[N][1]=j*a;
                   p[N][2]=k*a;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p[N][0], p[N][1], p[N][2]);
                   N=N+1;
                   p[N][0]=(2*i+1)*a/2;
                   p[N][1]=(2*j+1)*a/2;
                   p[N][2]=k*a;
                   N=N+1;
                   p[N][0]=i*a;
                   p[N][1]=(2*j+1)*a/2;
                   p[N][2]=(2*k+1)*a/2;
                   N=N+1;
                   p[N][0]=(2*i+1)*a/2;
                   p[N][1]=j*a;
                   p[N][2]=(2*k+1)*a/2;
                   N=N+1;
               }
           }
       }
    printf("Total number of atom N=%d\n", N);
    
    double va[3][3] = {{a*nx, 0, 0}, {0, a*ny, 0}, {0, 0, a*nz}};
    //calculation of vector b of given a defined by question
    i=0;
    j=1;
    k=2;
    for (int a=0; a<3; ++a) {
        cross_product(va[j][0], va[j][1], va[j][2], va[k][0], va[k][1], va[k][2]);
        for (int l=0; l<3; l++) {
            vb[i][l]=inter_cp[l];
            //printf("%lf\n", vb[i][l]);
        }
        dot_product(va[i][0], va[i][1], va[i][2], vb[i][0], vb[i][1], vb[i][2]);
        for (int l=0; l<3; l++) {
            vb[i][l]=vb[i][l]/inter_dp;
            //printf("%lf\n", vb[i][l]);
        }
        //printf("i=%d j=%d k=%d\n", i, j, k);
        h=k;
        k=i;
        i=j;
        j=h;
    }
    double *** vn = (double ***)malloc(N*sizeof(double**));
    for (i = 0; i< N; i++) {
        vn[i] = (double **) malloc(N*sizeof(double *));
            for (j = 0; j < N; j++) {
                vn[i][j] = (double *)malloc(3*sizeof(double));
            }
    }
    double  rn[N][3];   //rn here is just a array to store the vetor w.r.t. to the atom we choose
    for (k=0; k<N; k++) {
        for (j=0; j<N; j++) {
            for (int i=0; i<3; i++) {
                rn[j][i]=p[j][i]-p[k][i];
            }
            for (int i=0; i<3; i++) {
                dot_product(rn[j][0], rn[j][1], rn[j][2], vb[i][0], vb[i][1], vb[i][2]);
                vn[k][j][i]= inter_dp;
                //printf("%d %d %d %lf\n",k,j,i,vn[i][j][k]);
                if (vn[k][j][i]>0.5) {
                    vn[k][j][i]=vn[k][j][i]-1;
                }
                if (vn[k][j][i]<-0.5) {
                    vn[k][j][i]=vn[k][j][i]+1;
                }
                //printf("%d %d %d %lf\n",k,j,i,vn[i][j][k]);
            }
        }
    }
    
    double *** r = (double ***)malloc(N*sizeof(double**));
    for (i = 0; i< N; i++) {
        r[i] = (double **) malloc(N*sizeof(double *));
        for (j = 0; j < N; j++) {
            r[i][j] = (double *)malloc(3*sizeof(double));
        }
    }
    
    double **distance;
    distance = (double **) malloc (sizeof(double *)*NoA);
    for (i=0; i<NoA; i++) {
        distance[i] = (double *)malloc(sizeof(double)*NoA);
    }
    
    for (k=0; k<N; k++) {
        //printf("This is for the %d-th atom\n", k);
        for (j=k+1; j<N; j++) {
            for (int i=0; i<3; i++) {
            r[k][j][i]=vn[k][j][0]*va[0][i]+vn[k][j][1]*va[1][i]+vn[k][j][2]*va[2][i];
            }
            //printf("N=%d r1=%lf r2=%lf r3=%lf\n", j, r[k][j][0], r[k][j][1], r[k][j][2]);
            distance[k][j]=sqrt(r[k][j][0]*r[k][j][0]+ r[k][j][1]*r[k][j][1] +r[k][j][2]*r[k][j][2]);
            distance[j][k]=distance[k][j];
            /*if (distance[k][j]<(10*sigma) && distance[k][j]>0 ) {
                printf("N=%d r=%lf \n", j, distance[k][j]);
            }*/
        }
    }
    
    //energy of each atom(for checking)
    /*total_energy=0;
    for (k=0; k<NoA; k++) {
        for (j=0; j<NoA; j++) {
            if (distance[k][j]<(10*sigma) && distance[k][j]>0) {
                total_energy=total_energy+lennard_jones(distance[k][j]);
            }
        }
    printf("total energy of %d-th atom E=%lf eV\n", k, total_energy);
    total_energy=0;
    }*/
    
    //calculation of total energy of Aluminium(as its structure is fcc)
    total_energy=0;
    for (k=0; k<N; k++) {
        for (j=k+1; j<N; j++) {
            if (distance[k][j]<(10*sigma)) {
                total_energy=total_energy+lennard_jones(distance[k][j]);
            }
        }
    }
    printf("total energy E=%lfeV\n", total_energy);
    printf("total energy per atom=%lfeV\n", total_energy/NoA);
    
    return 0;
}

void cross_product(double x1,double x2,double x3,double y1,double y2,double y3){
    double cross_p[3] = {x2*y3-x3*y2, -(x1*y3-x3*y1), x1*y2-x2*y1};
    for (int i=0; i<3; i++) {
        inter_cp[i] = cross_p[i];
    }
}

void dot_product(double x1,double x2,double x3,double y1,double y2,double y3){
    inter_dp=x1*y1+x2*y2+x3*y3;
}

double lennard_jones(double r){
    double y, t;
    t=sigma/r;
    y=4*epsilon*(pow(t, 12)-pow(t, 6));
    return y;
}
