//
//  main.c
//  4_part4
//
//  Created by Leon Kwok on 7/12/2019.
//  Copyright © 2019 Leon Kwok. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
double vb[3][3], inter_cp[3];
double inter_dp, h;
double sigma, epsilon;
int main() {
    double a;   //lattice constant
    int i,j,k,NoA, N;
    int nx, ny, nz;
    double total_energy, old_total_energy, sigma_a, alpha, gamma;
    void cross_product(double,double,double,double,double,double);
    void dot_product(double,double,double,double,double,double);
    double CB(double, int, int);
    double grad_CB(double, double, int, int);
    
    FILE *diffm;
    diffm=fopen("m_part4.csv", "w");
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
    
    double p[2*NoA][3];
    double temp_p[2*NoA][3];
    for (i=0; i<NoA; i++) {
        p[2*i][0]=p_Mg[i][0];
        p[2*i][1]=p_Mg[i][1];
        p[2*i][2]=p_Mg[i][2];
        //printf("Na\t N=%d\t %lf\t %lf\t %lf\n",2*i, position[2*i][0], position[2*i][1], position[2*i][2]);
        p[2*i+1][0]=p_O[i][0];
        p[2*i+1][1]=p_O[i][1];
        p[2*i+1][2]=p_O[i][2];
        //printf("Cl\t N=%d\t %lf\t %lf\t %lf\n",2*i+1, position[2*i+1][0], position[2*i+1][1], position[2*i+1][2]);
    }
    N=2*N;
    //printf("N=%d\n", N);
    
    double min, max, ran_NO;
    min = -.05;
    max = .05;
    srand(756213);
    for (int i=0; i<N; i++) {
        //printf("before before x=%lf y=%lf z=%lf\n", p[i][0], p[i][1], p[i][2]);
        for (j=0; j<3; j++) {
            //printf("p1=%lf\n", p[i][j]);
            ran_NO=(max - min)*rand()/(RAND_MAX+1.0)+min;
            //printf("ran=%lf\n", ran_NO);
            p[i][j]=p[i][j]+ran_NO*a;
            //printf("p2=%lf\n", p[i][j]);
        }
        //printf("before x=%lf y=%lf z=%lf\n", p[i][0], p[i][1], p[i][2]);
    }
    
    double ***r = (double ***)malloc(N*sizeof(double**));
    for (i = 0; i< N; i++) {
        r[i] = (double **) malloc(N*sizeof(double *));
        for (j = 0; j < N; j++) {
            r[i][j] = (double *)malloc(3*sizeof(double));
        }
    }
    
    double **distance;
    distance = (double **) malloc (sizeof(double *)*2*NoA);
    for (i=0; i<2*NoA; i++) {
        distance[i] = (double *)malloc(sizeof(double)*2*NoA);
    }
    
    double *** vn = (double ***)malloc(N*sizeof(double**));
    for (i = 0; i< N; i++) {
        vn[i] = (double **) malloc(N*sizeof(double *));
            for (j = 0; j < N; j++) {
                vn[i][j] = (double *)malloc(3*sizeof(double));
            }
    }
    
    double ***r_f = (double ***)malloc(N*sizeof(double**));
    for (i = 0; i< N; i++) {
        r_f[i] = (double **) malloc(N*sizeof(double *));
        for (j = 0; j < N; j++) {
            r_f[i][j] = (double *)malloc(3*sizeof(double));
        }
    }
    
    double **distance_f;
    distance_f = (double **) malloc (sizeof(double *)*2*NoA);
    for (i=0; i<2*NoA; i++) {
        distance_f[i] = (double *)malloc(sizeof(double)*2*NoA);
    }
    
    double ***g1 = (double ***)malloc(N*sizeof(double**));
    for (i = 0; i< N; i++) {
        g1[i] = (double **) malloc(N*sizeof(double *));
        for (j = 0; j < N; j++) {
            g1[i][j] = (double *)malloc(3*sizeof(double));
        }
    }
    
    double ***g2 = (double ***)malloc(N*sizeof(double**));
    for (i = 0; i< N; i++) {
        g2[i] = (double **) malloc(N*sizeof(double *));
        for (j = 0; j < N; j++) {
            g2[i][j] = (double *)malloc(3*sizeof(double));
        }
    }
    
    double **g1_vector;
    g1_vector = (double **) malloc (sizeof(double *)*2*NoA);
    for (i=0; i<2*NoA; i++) {
        g1_vector[i] = (double *)malloc(sizeof(double)*3);
    }
    
    double **old_g1_vector;
    old_g1_vector = (double **) malloc (sizeof(double *)*2*NoA);
    for (i=0; i<2*NoA; i++) {
        old_g1_vector[i] = (double *)malloc(sizeof(double)*3);
    }
    
    double **g2_vector;
    g2_vector = (double **) malloc (sizeof(double *)*2*NoA);
    for (i=0; i<2*NoA; i++) {
        g2_vector[i] = (double *)malloc(sizeof(double)*3);
    }
    
    double **hh;
    hh = (double **) malloc (sizeof(double *)*2*NoA);
    for (i=0; i<2*NoA; i++) {
        hh[i] = (double *)malloc(sizeof(double)*3);
    }
    
    double va[3][3] = {{a*nx, 0, 0}, {0, a*ny, 0}, {0, 0, a*nz}};
    //calculation of vector b of given a defined by question
    i=0;
    j=1;
    k=2;
    for (int d=0; d<3; ++d) {
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
    
    double  rn[N][3];   //rn here is just a array to store the vetor w.r.t. to the atom we choose
    int m=0;
    while(m<100){
        printf("m=%d\n", m);
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
        
        for (k=0; k<N; k++) {
            //printf("This is for the %d-th atom\n", k);
            for (j=0; j<N; j++) {
                for (int i=0; i<3; i++) {
                r[k][j][i]=vn[k][j][0]*va[0][i]+vn[k][j][1]*va[1][i]+vn[k][j][2]*va[2][i];
                }
                distance[k][j]=sqrt(r[k][j][0]*r[k][j][0]+ r[k][j][1]*r[k][j][1] +r[k][j][2]*r[k][j][2]);
            }
        }
        
        old_total_energy=total_energy;
        total_energy=0;
        for (k=0; k<N; k++) {
            for (j=k+1; j<N; j++) {
                if (distance[k][j]<(20*a)) {
                    total_energy=total_energy+CB(distance[k][j], k, j);
                }
            }
        }
        printf("total energy E=%lfeV\n", total_energy);
        printf("total energy per atom=%lfeV\n", total_energy/(2*NoA));
        fprintf(diffm, "%d\t %lf\n", m, total_energy/(2*NoA));
        //printf("old=%lf\t new=%lf\n", old_total_energy, total_energy);
        //printf("per_change=%lf\n", (total_energy-old_total_energy)/old_total_energy*100);
        if (fabs((total_energy-old_total_energy)/old_total_energy*100)<1e-6) {
            //printf("total energy E=%lfeV\n", total_energy);
            //printf("total energy per atom=%lfeV\n", total_energy/NoA);
            break;
        }
        
        for (k=0; k<N; k++) {
            for (i=0; i<3; i++) {
                old_g1_vector[k][i]=g1_vector[k][i];
            }
        }
        
        for (k=0; k<N; k++) {
            for (j=0; j<N ; j++) {
                if (k==j) {
                    continue;
                }
                for (i=0; i<3; i++) {
                    g1[k][j][i] = -grad_CB(distance[k][j], -r[k][j][i], k, j);
                    //printf("%d %d g1=%lf\n",k,j, g1[k][j][i]);
                }
            }
        }
        
        for (k=0; k<N; k++) {
            for (i=0; i<3; i++) {
                g1_vector[k][i]=0;
            }
            for (j=0; j<N ; j++) {
                if (k==j) {
                    continue;
                }
                for (i=0; i<3; i++) {
                    g1_vector[k][i]+=g1[k][j][i];
                }
            }
        }
        
        sigma_a=1e-5;
        gamma=0;
        double scalar1=0,scalar2=0;
        for (k=0; k<N; k++) {
            if (m==0) {
                continue;
            }
            else{
                for (i=0; i<3; i++) {
                    scalar1-=g1_vector[k][i]*(g1_vector[k][i]-old_g1_vector[k][i]);
                    scalar2+=(old_g1_vector[k][i]*old_g1_vector[k][i]);
                }
                //gamma[k][i]=(g1_vector[k][i]*g1_vector[k][i])/(old_g1_vector[k][i]*old_g1_vector[k][i]);
            }
        }
        gamma=scalar1/scalar2;
        if (gamma>=0) {
            gamma=gamma;
        } else {
            gamma=0;
        }
        for (k=0; k<N; k++) {
            if (m==0) {
                for (i=0; i<3; i++) {
                    hh[k][i]=g1_vector[k][i];
                }
            } else {
                for (i=0; i<3; i++) {
                    hh[k][i]=g1_vector[k][i]+gamma*hh[k][i];
                }
            }
        }
        
        for (k=0; k<N; k++) {
            for (i=0; i<3; i++) {
                temp_p[k][i]=p[k][i]+sigma_a*hh[k][i];
            }
        }
        
        for (k=0; k<N; k++) {
            for (j=0; j<N; j++) {
                for (int i=0; i<3; i++) {
                    rn[j][i]=temp_p[j][i]-temp_p[k][i];
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
        
        for (k=0; k<N; k++) {
            //printf("This is for the %d-th atom\n", k);
            for (j=0; j<N; j++) {
                for (int i=0; i<3; i++) {
                r_f[k][j][i]=vn[k][j][0]*va[0][i]+vn[k][j][1]*va[1][i]+vn[k][j][2]*va[2][i];
                }
                distance_f[k][j]=sqrt(r_f[k][j][0]*r_f[k][j][0]+r_f[k][j][1]*r_f[k][j][1]+r_f[k][j][2]*r_f[k][j][2]);
            }
        }
        
        for (k=0; k<N; k++) {
            for (j=0; j<N; j++) {
                if (k==j) {
                    continue;
                }
                for (i=0; i<3; i++) {
                    g2[k][j][i] = grad_CB(distance_f[k][j], -r_f[k][j][i], k, j);
                    //printf("%d %d g2=%lf\n",k,j, g2[k][j][i]);
                    //printf("%d %d g2=%lf\n",k,j, g2[k][j][i]);
                }
            }
        }
        
        for (k=0; k<N; k++) {
             for (i=0; i<3; i++) {
                 g2_vector[k][i]=0;
             }
             for (j=0; j<N ; j++) {
                 if (k==j) {
                     continue;
                 }
                 for (i=0; i<3; i++) {
                     g2_vector[k][i]+=g2[k][j][i];
                 }
             }
         }
        
        alpha=0;
        scalar1=0;scalar2=0;
        for (k=0; k<N; k++) {
                for (i=0; i<3; i++) {
                    scalar1-=g1_vector[k][i]*hh[k][i];
                    scalar2+=(g2_vector[k][i]+g1_vector[k][i])*hh[k][i];
                    //printf("%f %f\n",scalar1,scalar2);
                }
        }
        alpha=-sigma_a*scalar1/scalar2;
        printf("alpha=%lf\n", alpha);
        for (k=0; k<N; k++) {
            for (i=0; i<3; i++) {
                p[k][i]=p[k][i]+alpha*hh[k][i];
            }
        }
        m++;
    }
    fclose(diffm);
    
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
    if(r>(10.0))return 0.0;
    y=A*exp(-r/B)-C*pow(r, -6)+(E*14.3981/r);
    //D*pow(10, -8)*pow(r, -8);
    //printf("%d %d %f, A=%f,\tC=%f,\t E=%f\n",k,j,r,A*exp(-B*r),-C*pow(10, -6)*pow(r, -6), (E*1.43981/r));
    return y;
}


double grad_CB(double r, double xyz, int k, int j){
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
    //if(r>11)return 0;
    if(r>(10.0))return 0.0;
    y=(-A/B*exp(-r/B)+6*C*pow(r, -7)-(E*14.3981/(r*r)))*xyz/r;
    //printf("%f %f %d %d\n",xyz,r,k,j);
    //D*pow(10, -8)*pow(r, -8);
    //printf("%d %d %f, A=%f,\tC=%f,\t E=%f\n",k,j,r,A*exp(-B*r),-C*pow(10, -6)*pow(r, -6), (E*1.43981/r));
    return y;
}
