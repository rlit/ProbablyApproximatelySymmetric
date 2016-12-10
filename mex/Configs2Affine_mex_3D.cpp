#include "mex.h"
#include <math.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    int numConfigs;
    
    // input variables
    double *configs;
    double tx,ty,tz,s,incl,theta,role;
    int sourceH,sourceW,sourceD;
    int targetH,targetW,targetD;
    int r1x, r1y,r1z, r2x, r2y, r2z;
    
    // output variables
    double *affines;
    int *insiders;
    
    
    /* Find the dimensions of the data */
    numConfigs = mxGetN(prhs[0]);
    
    
    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(12, numConfigs, mxREAL );
    plhs[1] = mxCreateNumericMatrix (1, numConfigs, mxINT32_CLASS, 0);
    
    
    /* Retrieve the input data */
    configs = mxGetPr(prhs[0]);
    sourceH = mxGetScalar(prhs[1]);
    sourceW = mxGetScalar(prhs[2]);
    sourceD = mxGetScalar(prhs[3]);
    targetH = mxGetScalar(prhs[4]);
    targetW = mxGetScalar(prhs[5]);
    targetD = mxGetScalar(prhs[6]);
    r1x = mxGetScalar(prhs[7]);
    r1y = mxGetScalar(prhs[8]);
    r1z = mxGetScalar(prhs[9]);
    r2x = mxGetScalar(prhs[10]);
    r2y = mxGetScalar(prhs[11]);
    r2z = mxGetScalar(prhs[12]);
    
    
    /* Create a pointer to the output data */
    affines = mxGetPr(plhs[0]);
    insiders = (int*)mxGetPr(plhs[1]);
    
    
    // MAIN LOOP
    for (int i = 0 ; i < numConfigs ; i++)
    {
        //if (i==48195)
        //	mexPrintf("MAIN LOOP\n");
// 		if (i%100000==0)
// 			mexPrintf("MAIN LOOP: config %d out of %d\n",i+1,numConfigs);
        
        tx = configs[7*i];
        ty = configs[7*i+1];
        tz = configs[7*i+2];
        s = configs[7*i+3];
        incl = configs[7*i+4];
        theta = configs[7*i+5];
        role = configs[7*i+6];
        
        // these are the coordinates of the unit sphere point, around which a rotation by 'role' is performed
        double ux = sin(incl)*cos(theta);
        double uy = sin(incl)*sin(theta);
        double uz = cos(incl);
        
        double sinRole = sin(role);
        double cosRole = cos(role);
        
//                 A = [cosRole+ux*ux*(1-cosRole),       ux*uy*(1-cosRole)-uz*sinRole,    ux*uz*(1-cosRole)+uy*sinRole,  tx;
//                         uy*ux*(1-cosRole)+uz*sinRole, cosRole+uy*uy*(1-cosRole),          uy*uz*(1-cosRole)-ux*sinRole,  ty;
//                         uz*ux*(1-cosRole)-uy*sinRole, uz*uy*(1-cosRole)+ux*sinRole,     cosRole+uz*uz*(1-cosRole)     ,  tz];
        
//                                                                         configs(gridInd,:) = [tx,ty,tz,s,incl,theta,role];
        
        if (s < 0)
        {
// the case of reflection w.r.t. the normal [ux,uy,uz]:
            // A =   [1-2uxux  -2uxuy  -2uxuz ;
            //       [-2uxuy   1-2uyuy -2uyuz ;
            //       [-2uxuz   -2uyuz  1-2uzuz]
            affines[12*i]   = 1 - 2*ux*ux;
            affines[12*i+1] = -2*ux*uy;
            affines[12*i+2] = -2*ux*uz;
            affines[12*i+3] = tx;
            affines[12*i+4] = -2*uy*ux;
            affines[12*i+5] = 1-2*uy*uy;
            affines[12*i+6] = -2*uy*uz;
            affines[12*i+7] = ty;
            affines[12*i+8] = -2*uz*ux;
            affines[12*i+9] = -2*uz*uy;
            affines[12*i+10] = 1-2*uz*uz;
            affines[12*i+11] = tz;
        }
        else
        {
            //syms tx ty r2 sx sy r1;
            // A =     [ sx*cos(r1)*cos(r2) - sy*sin(r1)*sin(r2), - sx*cos(r1)*sin(r2) - sy*cos(r2)*sin(r1), tx]
            //         [ sx*cos(r2)*sin(r1) + sy*cos(r1)*sin(r2),   sy*cos(r1)*cos(r2) - sx*sin(r1)*sin(r2), ty]
            //         [                                       0,                                         0,  1]
            //        matrixConfigs(i,:) = [A(1,1),A(1,2),A(1,3),A(2,1),A(2,2),A(2,3)];
            
            affines[12*i]   = s*(cosRole+ux*ux*(1-cosRole));
            affines[12*i+1] = s*(ux*uy*(1-cosRole)-uz*sinRole);
            affines[12*i+2] = s*(ux*uz*(1-cosRole)+uy*sinRole);
            affines[12*i+3] = tx;
            affines[12*i+4] = s*(uy*ux*(1-cosRole)+uz*sinRole);
            affines[12*i+5] = s*(cosRole+uy*uy*(1-cosRole));
            affines[12*i+6] = s*(uy*uz*(1-cosRole)-ux*sinRole);
            affines[12*i+7] = ty;
            affines[12*i+8] = s*(uz*ux*(1-cosRole)-uy*sinRole);
            affines[12*i+9] = s*(uz*uy*(1-cosRole)+ux*sinRole);
            affines[12*i+10] = s*(cosRole+uz*uz*(1-cosRole));
            affines[12*i+11] = tz;
        }
        //mexPrintf("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",a11,a12,a13,a21,a22,a23);
        
        // check if 8 corners of box are 'almost' within bounds
        double c1x = affines[12*i]  *(1-(r1x+1)) + affines[12*i+1]*(1-(r1y+1)) + affines[12*i+2]*(1-(r1z+1)) + (r2x+1) + tx;
        double c1y = affines[12*i+4]*(1-(r1x+1)) + affines[12*i+5]*(1-(r1y+1)) + affines[12*i+6]*(1-(r1z+1)) + (r2y+1) + ty;
        double c1z = affines[12*i+8]*(1-(r1x+1)) + affines[12*i+9]*(1-(r1y+1)) + affines[12*i+10]*(1-(r1z+1)) + (r2z+1) + tz;
        double c2x = affines[12*i]  *(sourceW-(r1x+1)) + affines[12*i+1]*(1-(r1y+1)) + affines[12*i+2]*(1-(r1z+1)) + (r2x+1) + tx;
        double c2y = affines[12*i+4]*(sourceW-(r1x+1)) + affines[12*i+5]*(1-(r1y+1)) + affines[12*i+6]*(1-(r1z+1)) + (r2y+1) + ty;
        double c2z = affines[12*i+8]*(sourceW-(r1x+1)) + affines[12*i+9]*(1-(r1y+1)) + affines[12*i+10]*(1-(r1z+1)) + (r2z+1) + tz;
        double c3x = affines[12*i]  *(sourceW-(r1x+1)) + affines[12*i+1]*(sourceH-(r1y+1)) + affines[12*i+2]*(1-(r1z+1)) + (r2x+1) + tx;
        double c3y = affines[12*i+4]*(sourceW-(r1x+1)) + affines[12*i+5]*(sourceH-(r1y+1)) + affines[12*i+6]*(1-(r1z+1)) + (r2y+1) + ty;
        double c3z = affines[12*i+8]*(sourceW-(r1x+1)) + affines[12*i+9]*(sourceH-(r1y+1)) + affines[12*i+10]*(1-(r1z+1)) + (r2z+1) + tz;
        double c4x = affines[12*i]  *(1-(r1x+1)) + affines[12*i+1]*(sourceH-(r1y+1)) + affines[12*i+2]*(1-(r1z+1)) + (r2x+1) + tx;
        double c4y = affines[12*i+4]*(1-(r1x+1)) + affines[12*i+5]*(sourceH-(r1y+1)) + affines[12*i+6]*(1-(r1z+1)) + (r2y+1) + ty;
        double c4z = affines[12*i+8]*(1-(r1x+1)) + affines[12*i+9]*(sourceH-(r1y+1)) + affines[12*i+10]*(1-(r1z+1)) + (r2z+1) + tz;
        
        double c5x = affines[12*i]  *(1-(r1x+1)) + affines[12*i+1]*(1-(r1y+1)) + affines[12*i+2]*(sourceD-(r1z+1)) + (r2x+1) + tx;
        double c5y = affines[12*i+4]*(1-(r1x+1)) + affines[12*i+5]*(1-(r1y+1)) + affines[12*i+6]*(sourceD-(r1z+1)) + (r2y+1) + ty;
        double c5z = affines[12*i+8]*(1-(r1x+1)) + affines[12*i+9]*(1-(r1y+1)) + affines[12*i+10]*(sourceD-(r1z+1)) + (r2z+1) + tz;
        double c6x = affines[12*i]  *(sourceW-(r1x+1)) + affines[12*i+1]*(1-(r1y+1)) + affines[12*i+2]*(sourceD-(r1z+1)) + (r2x+1) + tx;
        double c6y = affines[12*i+4]*(sourceW-(r1x+1)) + affines[12*i+5]*(1-(r1y+1)) + affines[12*i+6]*(sourceD-(r1z+1)) + (r2y+1) + ty;
        double c6z = affines[12*i+8]*(sourceW-(r1x+1)) + affines[12*i+9]*(1-(r1y+1)) + affines[12*i+10]*(sourceD-(r1z+1)) + (r2z+1) + tz;
        double c7x = affines[12*i]  *(sourceW-(r1x+1)) + affines[12*i+1]*(sourceH-(r1y+1)) + affines[12*i+2]*(sourceD-(r1z+1)) + (r2x+1) + tx;
        double c7y = affines[12*i+4]*(sourceW-(r1x+1)) + affines[12*i+5]*(sourceH-(r1y+1)) + affines[12*i+6]*(sourceD-(r1z+1)) + (r2y+1) + ty;
        double c7z = affines[12*i+8]*(sourceW-(r1x+1)) + affines[12*i+9]*(sourceH-(r1y+1)) + affines[12*i+10]*(sourceD-(r1z+1)) + (r2z+1) + tz;
        double c8x = affines[12*i]  *(1-(r1x+1)) + affines[12*i+1]*(sourceH-(r1y+1)) + affines[12*i+2]*(sourceD-(r1z+1)) + (r2x+1) + tx;
        double c8y = affines[12*i+4]*(1-(r1x+1)) + affines[12*i+5]*(sourceH-(r1y+1)) + affines[12*i+6]*(sourceD-(r1z+1)) + (r2y+1) + ty;
        double c8z = affines[12*i+8]*(1-(r1x+1)) + affines[12*i+9]*(sourceH-(r1y+1)) + affines[12*i+10]*(sourceD-(r1z+1)) + (r2z+1) + tz;
        
        // allow to exceed boundary by at most a factor of 5
//         int exceededX = sourceW/2;
//         int exceededY = sourceH/2;
//         int exceededZ = sourceD/2;
// 		int k = int((c1x>-exceededX)&&(c1x<targetW+exceededX)&&(c1y>-exceededY)&&(c1y<targetH+exceededY)&&(c1z>-exceededZ)&&(c1z<targetD+exceededZ)&&
//                     (c2x>-exceededX)&&(c2x<targetW+exceededX)&&(c2y>-exceededY)&&(c2y<targetH+exceededY)&&(c2z>-exceededZ)&&(c2z<targetD+exceededZ)&&
//                     (c3x>-exceededX)&&(c3x<targetW+exceededX)&&(c3y>-exceededY)&&(c3y<targetH+exceededY)&&(c3z>-exceededZ)&&(c3z<targetD+exceededZ)&&
//                     (c4x>-exceededX)&&(c4x<targetW+exceededX)&&(c4y>-exceededY)&&(c4y<targetH+exceededY)&&(c4z>-exceededZ)&&(c4z<targetD+exceededZ)&&
//                     (c5x>-exceededX)&&(c5x<targetW+exceededX)&&(c5y>-exceededY)&&(c5y<targetH+exceededY)&&(c5z>-exceededZ)&&(c5z<targetD+exceededZ)&&
//                     (c6x>-exceededX)&&(c6x<targetW+exceededX)&&(c6y>-exceededY)&&(c6y<targetH+exceededY)&&(c6z>-exceededZ)&&(c6z<targetD+exceededZ)&&
//                     (c7x>-exceededX)&&(c7x<targetW+exceededX)&&(c7y>-exceededY)&&(c7y<targetH+exceededY)&&(c7z>-exceededZ)&&(c7z<targetD+exceededZ)&&
//                     (c8x>-exceededX)&&(c8x<targetW+exceededX)&&(c8y>-exceededY)&&(c8y<targetH+exceededY)&&(c8z>-exceededZ)&&(c8z<targetD+exceededZ));
        //int kk = (5==5);
        //int kkk = (5==4);
        //if (i%100==0)
        //	mexPrintf("%d,%d,%d\n",k,kk,kkk);
// 		if(k)
// 			mexPrintf("good at config %d\n",i);
        int k = 1;
        insiders[i] = k;
    }
    
}
//
//%
//cornersX = [1 w1 w1 1];
//cornersY = [1 1 h1 h1];
//cornersA = round(a2x2*[cornersX-(r1x+1);cornersY-(r1y+1)]);
//cornerAxs = round(cornersA(1,:) + (r2x+1)  + a(1,3));
//cornerAys = round(cornersA(2,:) + (r2y+1)  + a(2,3));
