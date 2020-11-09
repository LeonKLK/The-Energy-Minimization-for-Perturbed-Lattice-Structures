//
//  main.c
//  3_part2
//
//  Created by Leon Kwok on 11/10/2019.
//  Copyright © 2019 Leon Kwok. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double vb[3][3], inter_cp[3];
double inter_dp, h;
//double sigma, epsilon;
int NoA;
int main() {
    double a;   //lattice constant
    int i,j,k,N;
    int nx, ny, nz;
    double total_energy;
    void cross_product(double,double,double,double,double,double);
    void dot_product(double,double,double,double,double,double);
    double CB(double, int, int);

    nx=5;
    ny=5;
    nz=5;
    a=4.212;    //using angstrom as unit
    NoA=nx*ny*nz*4;

    double p_Mg[NoA][3];
    N=0;
    for (k=0; k<nz; ++k) {
           for (j=0; j<ny; ++j) {
               for (i=0; i<nx; ++i) {
                   p_Mg[N][0]=i*a;
                   p_Mg[N][1]=j*a;
                   p_Mg[N][2]=k*a;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_Mg[N][0], p_Mg[N][1], p_Mg[N][2]);
                   N=N+1;
                   p_Mg[N][0]=(2*i+1)*a/2;
                   p_Mg[N][1]=(2*j+1)*a/2;
                   p_Mg[N][2]=k*a;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_Mg[N][0], p_Mg[N][1], p_Mg[N][2]);
                   N=N+1;
                   p_Mg[N][0]=i*a;
                   p_Mg[N][1]=(2*j+1)*a/2;
                   p_Mg[N][2]=(2*k+1)*a/2;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_Mg[N][0], p_Mg[N][1], p_Mg[N][2]);
                   N=N+1;
                   p_Mg[N][0]=(2*i+1)*a/2;
                   p_Mg[N][1]=j*a;
                   p_Mg[N][2]=(2*k+1)*a/2;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_Mg[N][0], p_Mg[N][1], p_Mg[N][2]);
                   N=N+1;
               }
           }
       }
    double p_O[NoA][3];
    N=0;
    for (k=0; k<nz; ++k) {
           for (j=0; j<ny; ++j) {
               for (i=0; i<nx; ++i) {
                   p_O[N][0]=i*a+a/2;
                   p_O[N][1]=j*a;
                   p_O[N][2]=k*a;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_O[N][0], p_O[N][1], p_O[N][2]);
                   N=N+1;
                   p_O[N][0]=(2*i+1)*a/2+a/2;
                   p_O[N][1]=(2*j+1)*a/2;
                   p_O[N][2]=k*a;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_O[N][0], p_O[N][1], p_O[N][2]);
                   N=N+1;
                   p_O[N][0]=i*a+a/2;
                   p_O[N][1]=(2*j+1)*a/2;
                   p_O[N][2]=(2*k+1)*a/2;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_O[N][0], p_O[N][1], p_O[N][2]);
                   N=N+1;
                   p_O[N][0]=(2*i+1)*a/2+a/2;
                   p_O[N][1]=j*a;
                   p_O[N][2]=(2*k+1)*a/2;
                   //printf("N=%d\t %lf\t %lf\t %lf\n",N, p_O[N][0], p_O[N][1], p_O[N][2]);
                   N=N+1;
               }
           }
       }
    /*N=0;
    for (i=0; i<NoA; ++i) {
        printf("N=%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",N, p_Mg[i][0], p_Mg[i][1], p_Mg[i][2], p_O[i][0], p_O[i][1], p_O[i][2]);
        N=N+1;
    }*/
    double position[2*NoA][3];
    for (i=0; i<NoA; i++) {
        position[2*i][0]=p_Mg[i][0];
        position[2*i][1]=p_Mg[i][1];
        position[2*i][2]=p_Mg[i][2];
        //printf("Na\t N=%d\t %lf\t %lf\t %lf\n",2*i, position[2*i][0], position[2*i][1], position[2*i][2]);
        position[2*i+1][0]=p_O[i][0];
        position[2*i+1][1]=p_O[i][1];
        position[2*i+1][2]=p_O[i][2];
        //printf("Cl\t N=%d\t %lf\t %lf\t %lf\n",2*i+1, position[2*i+1][0], position[2*i+1][1], position[2*i+1][2]);
    }
    
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
    double *** vn = (double ***)malloc(2*NoA*sizeof(double**));
    for (i = 0; i<2*NoA; i++) {
        vn[i] = (double **) malloc(2*NoA*sizeof(double *));
            for (j = 0; j <2*NoA; j++) {
                vn[i][j] = (double *)malloc(3*sizeof(double));
            }
    }
    double  rn[2*NoA][3];   //rn here is just a array to store the vetor w.r.t. to the atom we choose
    for (k=0; k<2*NoA; k++) {
        for (j=0; j<2*NoA; j++) {
            for (int i=0; i<3; i++) {
                rn[j][i]=position[j][i]-position[k][i];
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
    
    double *** r = (double ***)malloc(2*NoA*sizeof(double**));
    for (i = 0; i<2*NoA; i++) {
        r[i] = (double **) malloc(2*NoA*sizeof(double *));
        for (j = 0; j <2*NoA; j++) {
            r[i][j] = (double *)malloc(3*sizeof(double));
        }
    }
    double **distance;
    distance = (double **) malloc (sizeof(double *)*2*NoA);
    for (i=0; i<2*NoA; i++) {
        distance[i] = (double *)malloc(sizeof(double)*2*NoA);
    }
    for (k=0; k<2*NoA; k++) {
        //printf("This is for the %d-th atom\n", k);
        for (j=k+1; j<2*NoA; j++) {
            for (int i=0; i<3; i++) {
            r[k][j][i]=vn[k][j][0]*va[0][i]+vn[k][j][1]*va[1][i]+vn[k][j][2]*va[2][i];
            }
            //printf("N=%d r1=%lf r2=%lf r3=%lf\n", j, r[k][j][0], r[k][j][1], r[k][j][2]);
            distance[k][j]=sqrt(r[k][j][0]*r[k][j][0]+ r[k][j][1]*r[k][j][1] +r[k][j][2]*r[k][j][2]);
            distance[j][k]=distance[k][j];
            /*if (distance[k][j]<(2) && distance[k][j]>0 ) {
                printf("N=%d r=%lf \n", j, distance[k][j]);
            }*/
        }
    }
    double cutoff_d=20*a;
    //energy of each atom
    /*total_energy=0;
    for (k=0; k<2*NoA; k++) {
        for (j=0; j<2*NoA; j++) {
            if (distance[k][j]<cutoff_d && distance[k][j]>0) {
                total_energy=total_energy+CB(distance[k][j], k, j);
            }
        }
    printf("total energy of %d-th atom has\t E=%lf eV\n", k, total_energy);
    total_energy=0;
    }*/
    
    //calculation of total energy of NaCl(as its structure is fcc)
    total_energy=0;
    for (k=0; k<2*NoA; k++) {
        for (j=k+1; j<2*NoA; j++) {
            if (distance[k][j]<(cutoff_d)) {
                total_energy=total_energy+CB(distance[k][j], k, j);
            }
        }
    }
    printf("total energy E=%lfeV\n", total_energy);
    printf("total energy per atom=%lfeV\n", total_energy/(2*NoA));
    
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

double CB(double r, int k, int j){
    double y, A, B, C, E;
    if ((k%2==0) && (j%2==0)) {     //Mg-Mg, even number
        A=0;
        B=0.4208;
        C=0;
        E=4;
    }
    else if ((k%2==1) && (j%2==1)) {     //O-O
        A=9547.96;
        B=0.21916;
        C=32;
        E=4;
    }
    else {                          //Mg-O
        A=1279.69;
        B=0.29969;
        C=0;
        E=-4;
    }
    
    //*8.9876e18/r)/1.6022e-19)
    y=A*exp(-r/B)-C*pow(r, -6)+(E*14.3981/r);
    //D*pow(10, -8)*pow(r, -8);
    //printf("%d %d %f, A=%f,\tC=%f,\t E=%f\n",k,j,r,A*exp(-B*r),-C*pow(10, -6)*pow(r, -6), (E*1.43981/r));
    return y;
}
