//
//  main.c
//  4_part1
//
//  Created by Leon Kwok on 28/11/2019.
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
    double total_energy, old_total_energy, sigma_a, alpha;
    void cross_product(double,double,double,double,double,double);
    void dot_product(double,double,double,double,double,double);
    double lennard_jones(double);
    double grad_LJ(double, double);
    
    sigma=2.6175;   //using angstroms as unit
    epsilon=4577.6*8.617333262145e-5;     //4577.6*boltzmann constant eV/K
    FILE *diffm;
    diffm=fopen("m_part1.csv", "w");
    nx=5;
    ny=5;
    nz=5;
    a=4.046;    //using angstroms as unit
    
    NoA=nx*ny*nz*4;
    
    double p[NoA][3];
    double temp_p[NoA][3];
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
                   p[N][2]=k*a*1.03;
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
    distance = (double **) malloc (sizeof(double *)*NoA);
    for (i=0; i<NoA; i++) {
        distance[i] = (double *)malloc(sizeof(double)*NoA);
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
    distance_f = (double **) malloc (sizeof(double *)*NoA);
    for (i=0; i<NoA; i++) {
        distance_f[i] = (double *)malloc(sizeof(double)*NoA);
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
    g1_vector = (double **) malloc (sizeof(double *)*NoA);
    for (i=0; i<NoA; i++) {
        g1_vector[i] = (double *)malloc(sizeof(double)*3);
    }
    
    double **g2_vector;
    g2_vector = (double **) malloc (sizeof(double *)*NoA);
    for (i=0; i<NoA; i++) {
        g2_vector[i] = (double *)malloc(sizeof(double)*3);
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
                if (distance[k][j]<(10*sigma)) {
                    total_energy=total_energy+lennard_jones(distance[k][j]);
                }
            }
        }
        printf("total energy E=%lfeV\n", total_energy);
        printf("total energy per atom=%lfeV\n", total_energy/NoA);
        fprintf(diffm, "%d\t %lf\n", m, total_energy/NoA);
        //printf("old=%lf\t new=%lf\n", old_total_energy, total_energy);
        //printf("per_change=%lf\n", (total_energy-old_total_energy)/old_total_energy*100);
        if (fabs((total_energy-old_total_energy)/old_total_energy*100)<0.005) {
            break;
        }
        
        for (k=0; k<N; k++) {
            for (j=0; j<N ; j++) {
                if (k==j) {
                    continue;
                }
                for (i=0; i<3; i++) {
                    g1[k][j][i] =-grad_LJ(distance[k][j], -r[k][j][i]);
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
            //printf("%d %f %f %f \n",k,g1_vector[k][0],g1_vector[k][1],g1_vector[k][2]);
        }
        
        sigma_a=1e-5;
        for (k=0; k<N; k++) {
            for (i=0; i<3; i++) {
                temp_p[k][i]=p[k][i]+sigma_a*g1_vector[k][i];
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
                    g2[k][j][i] = grad_LJ(distance_f[k][j], -r_f[k][j][i]);
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
        double scalar1=0,scalar2=0;
        for (k=0; k<N; k++) {
                for (i=0; i<3; i++) {
                    scalar1-=g1_vector[k][i]*g1_vector[k][i];
                    scalar2+=(g2_vector[k][i]+g1_vector[k][i])*g1_vector[k][i];
                    //printf("%f %f\n",scalar1,scalar2);
                }
            //printf("alpha=%lf\n", alpha);
        }
        alpha=-sigma_a*scalar1/scalar2;
        printf("alpha=%lf\n", alpha);
        for (k=0; k<N; k++) {
            for (i=0; i<3; i++) {
                p[k][i]=p[k][i]+alpha*g1_vector[k][i];
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

double lennard_jones(double r){
    double y, t;
    t=sigma/r;
    y=4*epsilon*(pow(t, 12)-pow(t, 6));
    return y;
}

double grad_LJ(double r, double xyz){
    double y, t;
    t=sigma/r;
    y=4*epsilon*(6*pow(t, 7)/r-12*pow(t, 13)/r)*xyz/r;
    return y;
}
