#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "omp.h"

#define nproc 4
#define FLOAT double
void update_eigensys(const FLOAT *pos, int len, FLOAT* axes, FLOAT *vecs);
int *sort_evals(gsl_vector* eval);
void get_shape(const void * posv, int len, void * vecsv, void * axesv);

/*
 * Calculates the shape of the halo via the inertia tensor according to Allgood's method
 * In      : posv  = The values of position of DM particles (const)
 * In      : len   = The number of DM particles
 * Modifies: vecsv = [v1,v2,v3] the eigenvectors of the shape vector
 * Modifies: axes  = [a, b, c ] the eigenvalues of the shape vector
 * Initial values: a,b,c = rad_________vecs = canonical base
 * 
 */
void get_shape(const void * posv, int len, void * vecsv, void * axesv) {
    
    omp_set_num_threads(nproc);
    // Indices
    int i = 0;

    // Interprets data 
    const FLOAT * pos = (FLOAT *) posv;
    FLOAT * vecs = (FLOAT *) vecsv;
    FLOAT * axes = (FLOAT *) axesv;
    
    /*
    printf("vec1 = %f  %f  %f\n", vecs[0],vecs[1],vecs[2]);
    printf("vec2 = %f  %f  %f\n", vecs[3],vecs[4],vecs[5]);
    printf("vec3 = %f  %f  %f\n", vecs[6],vecs[7],vecs[8]);
    printf("axs = %f  %f  %f\n", axes[0],axes[1],axes[2]);
    */

    // Keeps eigenvalues for convergence calculations
    FLOAT a = axes[0];
    FLOAT b = axes[0];
    FLOAT c = axes[0];
    FLOAT conv = 100;
    do{
        i++;
	//printf("Iteration %d ______________________________________________________________\n", i);
        update_eigensys(pos,len,axes,vecs);
	conv = (1./3.)*((fabs(a-axes[0])/a)+(fabs(b-axes[1])/b)+(fabs(c-axes[2])/c));
	//printf("Axes = %f  %f  %f\n", a,b,c);
	//printf("Convergence = %f  %f\n",conv,pow(10.0,-6));
	a = axes[0];
	b = axes[1];
	c = axes[2];
	if( a < 0 ){
	    break;
	}
    }while( fabs(conv) >= pow(10.0,-6) );
    //printf("Convergence achived at:   %d\n",i);
}

/* 
 * Diagonalizes the inertia tensor and update eigen vals and eigenvecs
 * vecs define the major axis of the ellipsoid
 * axes define the axes sizes
 *
 */
void update_eigensys(const FLOAT *pos, int len, FLOAT* axes, FLOAT *vecs){
    
    int i;
    int j;
    // Allocates the inertia matrix
    double *inertia = malloc(9*sizeof(double)); 

 
    for( i = 0; i < 9; i++ ){
        inertia[i] = 0.0;
    }

    // Number of particles enclosed by ellipse and temporal variables for inertia
    int nump = 0;
    FLOAT i1,i2,i3,i4,i5,i6;
    i1=0;i2=0;i3=0;i4=0;i5=0;i6=0;

    // Filters positions and calculates the inertia tensor
//#pragma omp parallel
{
    #pragma omp parallel for shared(pos,vecs,axes,len) reduction(+:nump,i1,i2,i3,i4,i5,i6)
    for( i = 0; i < len; i++ ){
 
        FLOAT val = pow((pos[3*i]*vecs[0]+pos[3*i+1]*vecs[1]+pos[3*i+2]*vecs[2])/axes[0],2);
	val      += pow((pos[3*i]*vecs[3]+pos[3*i+1]*vecs[4]+pos[3*i+2]*vecs[5])/axes[1],2);
	val      += pow((pos[3*i]*vecs[6]+pos[3*i+1]*vecs[7]+pos[3*i+2]*vecs[8])/axes[2],2);
	
	if( (val <= 1) && (val > 0) ){
	    nump += 1;
	    
	    // Inertia term for this particle is calculated triangular superior
	    //inertia[3*l+n] += (double) pos[3*i+l]*pos[3*i+n]/(pow(axes[0],2)*val);
	    i1 += (double) pos[3*i+0]*pos[3*i+0]/(pow(axes[0],2)*val);
	    i2 += (double) pos[3*i+0]*pos[3*i+1]/(pow(axes[0],2)*val);
	    i3 += (double) pos[3*i+0]*pos[3*i+2]/(pow(axes[0],2)*val);
	    i4 += (double) pos[3*i+1]*pos[3*i+1]/(pow(axes[0],2)*val);
	    i5 += (double) pos[3*i+1]*pos[3*i+2]/(pow(axes[0],2)*val);
	    i6 += (double) pos[3*i+2]*pos[3*i+2]/(pow(axes[0],2)*val);
	}		
    }
}
    // Actualizes inertia
    inertia[0] = i1;
    inertia[1] = i2;
    inertia[2] = i3;
    inertia[4] = i4;
    inertia[5] = i5;
    inertia[8] = i6;
        
    /*#pragma omp critical (INSIDE PRAGMA PARALLEL)
    int l;
    int n;
    for( l = 0; l < 3; l++ ){
        for( n = l; n < 3; n++ ){
	  inertia[3*l+n] = temp[3*l+n] 
	}
    }
    */

    // If there are less than 3000 particles, break
    if( nump < 3000){
        axes[0] = 0;
        return; 
    }

    // Inertia matrix is symmetrized
    int l;
    int n;
    
    // Calculates the inertia tensor for this particle
    for( l = 0; l < 3; l++ ){
	for( n = l; n < 3; n++){
	    inertia[3*n+l] = inertia[3*l+n];
	}
    }

    /*
    printf("\n_________________________INERTIA_____________________________\n\n");
    printf("vec1 = %f  %f  %f\n", inertia[0],inertia[1],inertia[2]);
    printf("vec1 = %f  %f  %f\n", inertia[3],inertia[4],inertia[5]);
    printf("vec1 = %f  %f  %f\n", inertia[6],inertia[7],inertia[8]);
    printf("\n---------------------------------------------------------------\n");
    */

    // Diagonalizes (symmetric) matrix GNU scientific
    gsl_matrix_view m = gsl_matrix_view_array(inertia, 3, 3);
    
    // Allocation
    gsl_vector *eval = gsl_vector_alloc(3);
    gsl_matrix *evec = gsl_matrix_alloc(3, 3);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
    // Eigensystem solve
    gsl_eigen_symmv(&m.matrix, eval, evec, w);
    // Sorts eigenvectors and eigenvalues
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);        
    // Free Workspace
    gsl_eigen_symmv_free(w);
    
    /*
    for (i = 0; i < 3; i++)
      {
        double eval_i  = gsl_vector_get(eval, i);
        //gsl_vector_view evec_i = gsl_matrix_column (evec, i);

        printf ("eigenvalue = %.4e\n", eval_i);
        printf ("eigenvector = \n");
	printf("[%f  %f  %f]\n",vecs[3*i],vecs[3*i+1],vecs[3*i+2]);
	  //gsl_vector_fprintf (stdout, 
	  //                &evec_i.vector, "%.4e");
      }
    */

    // Actualizes the eigensystem conserving biggest semiaxis
    FLOAT v0 = ((FLOAT) gsl_vector_get (eval, 0));
    for( i = 0; i < 3; i++ ){
        if( i != 0 ){
	    axes[i] = axes[0]*sqrt(fabs(((FLOAT) gsl_vector_get(eval,i))/v0));
	}
	for( j = 0; j < 3; j++ ){
	    vecs[3*i + j] = gsl_matrix_get(evec, j, i);
	}
    }
    
    // Free memory
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    free(inertia);
}

int main(){
  return 0;
}
