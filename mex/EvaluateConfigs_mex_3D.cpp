#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{
	// parameters
	int h1, w1, d1, h2, w2, d2;
	int r1x, r1y, r1z, r2x, r2y, r2z;

	int numPoints, numConfigs;
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;

	// input variables
	double *img1, *img2;
	double *affines; // the matrices of the transformations
	int *xs, *ys, *zs;
	
    // Should we use photometrics?
    int usePhotometrics;
        
	// helper variables
	int *xs_centered, *ys_centered, *zs_centered;
	double *valsI1;
	int targetPoint_x, targetPoint_y, targetPoint_z;
	int targetInd;
	//double score;
	int maxInd2img2;

	// output variables
	double *distances;


	//// checking inputs
	//if (nrhs != 5)
	//	mexErrMsgTxt("The number of input arguments must be 5.");
	//if (nlhs != 1)
	//	mexErrMsgTxt("The number of output arguments must be 5.");


	// INPUTS:           dimensions:                type:
	// ------	         -----------					----
	// 0] image I1       h1 x w1 x d1				double
	// 1] image I2       h2 x w2 x d2				double
	// 2] configs        numConfigs x 6				double
	// 3] xs             1 x numPoints				int32
	// 4] ys		     1 x numPoints				int32
	// 5] zs		     1 x numPoints				int32
	// 6] photoInv		 1 x 1						int32

	// OUTPUTS:		 dimensions:				type:
	// -------		 -----------					-----
	// distances	 1 x numConfigs				double


	/* Find the dimensions of the data */
	h1 = (mxGetDimensions(prhs[0]))[0];
	w1 = (mxGetDimensions(prhs[0]))[1];
	d1 = (mxGetDimensions(prhs[0]))[2];
	h2 = (mxGetDimensions(prhs[1]))[0];
	w2 = (mxGetDimensions(prhs[1]))[1];
	d2 = (mxGetDimensions(prhs[1]))[2];
// 	w1 = mxGetM(prhs[0]);
// 	h2 = mxGetN(prhs[1]);
// 	w2 = mxGetM(prhs[1]);
	numConfigs = mxGetN(prhs[2]);
	numPoints = mxGetN(prhs[3]);

	r1x = 0.5*(w1-1);
	r1y = 0.5*(h1-1);
	r1z = 0.5*(d1-1);
	r2x = 0.5*(w2-1);
	r2y = 0.5*(h2-1);
	r2z = 0.5*(d2-1);

	maxInd2img2 = d2*h2*w2 - 1;


	/* Create an mxArray for the output data */
	plhs[0] = mxCreateDoubleMatrix(1, numConfigs, mxREAL );

	/* Create an mxArrays for temporary data */
	xs_centered = (int *)malloc(numPoints*sizeof(int));
	ys_centered = (int *)malloc(numPoints*sizeof(int));
	zs_centered = (int *)malloc(numPoints*sizeof(int));
	valsI1 = (double *)malloc(numPoints*sizeof(double));
        
    /* Store target pixel locations - x and y */
    //double* xs_target = (double *)malloc(numPoints*sizeof(double));
	//double* ys_target = (double *)malloc(numPoints*sizeof(double));
	//double* zs_target = (double *)malloc(numPoints*sizeof(double));
            

	/* Retrieve the input data */
	img1 = mxGetPr(prhs[0]);
	double* tmp_img2 = mxGetPr(prhs[1]);
	affines = mxGetPr(prhs[2]);
	xs = (int*)mxGetPr(prhs[3]);
	ys = (int*)mxGetPr(prhs[4]);
	zs = (int*)mxGetPr(prhs[5]);
    usePhotometrics = mxGetScalar(prhs[6]);

	// img2 is of height 3*img2 (this padding is for not needing to check bounds)
	img2 = (double*)malloc(5*d2*h2*w2*sizeof(double));
	memset(img2,2,5*d2*h2*w2*sizeof(double));
	memcpy(img2+2*d2*h2*w2,tmp_img2,d2*h2*w2*sizeof(double));
	//

	int cN = mxGetN(prhs[2]);
	int cM = mxGetM(prhs[2]);


	/*Centered pointes*/
	for (int i = 0 ; i < numPoints ; i++)
	{
		xs_centered[i] = xs[i]-(r1x+1);
		ys_centered[i] = ys[i]-(r1y+1);
		zs_centered[i] = zs[i]-(r1z+1);
	}

	/*Precalculating source point indices into I1 (and the values themselves)*/
	for (int j = 0; j < numPoints ; j++)
	{
		//sourceInds[j] = (ys[j] - 1)*w1 + xs[j];
		valsI1[j] = img1[(zs[j] - 1)*h1*w1+ (xs[j] - 1)*h1 + ys[j]-1]; // -1 is for c
	}


	/* Create a pointer to the output data */
	distances = mxGetPr(plhs[0]);


	// MAIN LOOP
	double tx, ty, tz, score, tmp1, tmp2, tmp3, *ptrVals;		
	int i,j, *ptrXsc, *ptrYsc, *ptrZsc;
#pragma omp parallel for if (numConfigs > 1000) num_threads(omp_get_num_procs()) default(shared) private(j,i,score, tx,ty,tz, tmp1,tmp2,tmp3, ptrVals,ptrXsc,ptrYsc,ptrZsc, targetPoint_x,targetPoint_y,targetPoint_z,targetInd, a11,a12,a13,a21,a22,a23,a31,a32,a33)
	for (i = 0 ; i < numConfigs ; i++)
	{

		/*if (i%100000==0)
			mexPrintf("MAIN LOOP: config %d out of %d\n",i+1,numConfigs);
		if (i==120804)
			mexPrintf("MAIN LOOP: config %d out of %d\n",i+1,numConfigs);*/

		a11 = affines[12*i];
		a12 = affines[12*i+1];
		a13 = affines[12*i+2];
		a21 = affines[12*i+4];
		a22 = affines[12*i+5];
		a23 = affines[12*i+6];
		a31 = affines[12*i+8];
		a32 = affines[12*i+9];
		a33 = affines[12*i+10];

		tx = affines[12*i+3];
		ty = affines[12*i+7];
		tz = affines[12*i+11];

		score = 0;
	
		tmp1 = (r2x+1) + tx + 0.5;
		tmp2 = (r2y+1) + ty + 0.5;
		tmp3 = (r2z+1) + tz + 0.5 + 2*d2; // + 2*d2 is for 'padding' of img2
		ptrVals = valsI1;
		ptrXsc = xs_centered;
		ptrYsc = ys_centered;
		ptrZsc = zs_centered;
		for (j = 0; j < numPoints ; j++)
			{
				targetPoint_x = int(a11*(*ptrXsc) + a12*(*ptrYsc) + a13*(*ptrZsc) + tmp1); // includes rounding
				targetPoint_y = int(a21*(*ptrXsc) + a22*(*ptrYsc) + a23*(*ptrZsc) + tmp2); // includes rounding
				targetPoint_z = int(a31*(*ptrXsc) + a32*(*ptrYsc) + a33*(*ptrZsc) + tmp3); // includes rounding
                                
				// FIX THE TARGET IND!!!!!!!!!!!!!!!!!!!!!!! 
				targetInd = (targetPoint_z - 1)*h2*w2 + (targetPoint_x - 1)*h2 + targetPoint_y - 1; // -1 is for c
				score += (fabs((*ptrVals) - img2[targetInd])); 
                                
				ptrVals++;
				ptrXsc++;
				ptrYsc++;
				ptrZsc++;
			}
		distances[i] = score/numPoints;
	}		

	/* Free the allocated arrays */
	free(xs_centered);
	free(ys_centered);
	free(zs_centered);
	free(valsI1);
	free(img2);

}