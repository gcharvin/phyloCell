// Copyright 2006-2007 The MathWorks, Inc.
#include "mex.h"
#include "FifoPriorityQueue.h"
#include "matrix.h"
#include <math.h>

#define getMax(x,y)     ((x)>(y)?(x):(y))
/*
 * mwSize= int or long int
 * L=output of the functon (double)
 * I=imput of the funtion
 * M=marqueur (double)
 *_T=the tipe of priority, usualy , double
 */

struct ellipseVar{
    double MajorAxisLength;
    double MinorAxisLength;
    double Eccentricity;
    double Orientation;
    double Centroid[2];
    double f1[2];
    double f2[2];
};

template<typename _T>
        void compute_watershed(_T *I, double *M, mwSize N, double *L, const mwSize *size_I) {
    
    FifoPriorityQueue<mwSize, _T> queue(FifoPriorityItemCompareFcn<mwSize, _T>::LowestPriorityFirst);
    /*=================================================================
     *declaration initialisation */
    mwSize s1= size_I[0];
    mwSize s2= size_I[1];
    mwSize neigh[8];
    mwSize edge;
    mwSize p;
    mwSize count=0;
    bool haveNeigh=false;
    double *penalized; /*the penalized pixels*/
    double *deleted; /*the pixels reanalized (deleted)*/
    double *F; /*already treted pixels (L without queue)*/
    mwSize *countPixelsObj; /*the nombre of pixels in each region*/
    mwSize numObj=0; /*the number of regions*/
    mwSize obj=0;/*actual object*/
    mwSize *x;
    mwSize *y;/*x and y of the pixels of actual object*/
    
    /*ellipse variable declaration*/
    struct ellipseVar stats;
    double xbar;
    double ybar;
    double difx;
    double dify;
    double uxx; 
    double uyy; 
    double uxy;
    double common;
    double num;
    double den;
    const double pi=3.14159265;
    double A;
    double a;
    double alpha;
    
    
    for (mwSize p = 0; p < N; p++) {
        L[p] = M[p];
        numObj=getMax(numObj,L[p]);
        if (M[p] != 0) {
            queue.push(p, I[p]); /*add marqueurs to queue*/
        }
    }
    
    countPixelsObj=(mwSize *) mxCalloc(numObj, sizeof(mwSize));
    
    /*====================================================================
     *main WHILE */
    while (! queue.isEmpty() ) {
        //count++;
        p = queue.topData();
        _T  v = queue.topPriority();
        queue.pop();

        edge=p%s1;
        if ((p>=s1) && (p<=N-s1) && ((edge!=0) && (edge!=s1-1))) {
            neigh[0]=p-s1-1;
            neigh[1]=p-s1;
            neigh[2]=p-s1+1;
            neigh[3]=p-1;
            neigh[4]=p+1;
            neigh[5]=p+s1-1;
            neigh[6]=p+s1;
            neigh[7]=p+s1+1;
            haveNeigh=true;
        }
        
        if (haveNeigh){
            for (mwSize j = 0; j < 8; j++) {
                if (L[neigh[j]]==0) {
                    L[neigh[j]]=L[p];
                    queue.push(neigh[j], I[neigh[j]]);
                    countPixelsObj[obj]++;/*count all the pixel of that region*/
                }
            }
        }
        
        
    }
    
}


//////////////////////////////////////////////////////////////////////////////
template<typename _T>
        inline void do_nan_check(_T *F, mwSize num_elements) {
    for (mwSize p = 0; p < num_elements; p++) {
        if (mxIsNaN(F[p])) {
            mexErrMsgIdAndTxt("Images:watershed:expectedNonNaN",
                    "%s",
                    "Input image may not contain NaNs.");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

#define SIZE_MISMATCH_ID  "Images:watershed_meyer:sizeMismatch"
#define SIZE_MISMATCH_MSG "I must be the same size as L."
//////////////////////////////////////////////////////////////////////////////
// Do some input validation.  Validation of the type of I, as well as
// the connectivity argument, is done in toolbox/images/images/watershed.m.
//////////////////////////////////////////////////////////////////////////////
void check_inputs(int nrhs, const mxArray *prhs[]) {
    mwSize num_dims;
    
    /*
     * if (nrhs != 3)
     * {
     * mexErrMsgIdAndTxt("Images:watershed_meyer:invalidNumInputs",
     * "%s",
     * "WATERSHED_MEYER needs 3 input arguments.");
     * }
     */
    
    if (!mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt("Images:watershed_meyer:nondoubleMinimaLabels",
                "%s",
                "Third input argument to WATERSHED_MEYER must be double.");
    }
    
    num_dims = mxGetNumberOfDimensions(prhs[0]);
    if (num_dims != mxGetNumberOfDimensions(prhs[1])) {
        mexErrMsgIdAndTxt(SIZE_MISMATCH_ID, "%s", SIZE_MISMATCH_MSG);
    }
    
    const mwSize *size_I = mxGetDimensions(prhs[0]);
    const mwSize *size_L = mxGetDimensions(prhs[1]);
    for (mwSize k = 0; k < num_dims; k++) {
        
        if (size_I[k] != size_L[k]) {
            mexErrMsgIdAndTxt(SIZE_MISMATCH_ID, "%s", SIZE_MISMATCH_MSG);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////






void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    void *I;
    double *L;
    mwSize num_elements;
    const mwSize *input_size;
    mwSize ndims;
    mxClassID class_id;
    
    
    (void) nlhs;  /* unused parameter */
    
    check_inputs(nrhs, prhs);
    
    if (mxIsLogical(prhs[0])) {
        I = mxGetLogicals(prhs[0]);
    }
    else {
        I = mxGetData(prhs[0]);
    }
    num_elements = mxGetNumberOfElements(prhs[0]);
    class_id = mxGetClassID(prhs[0]);
    input_size = mxGetDimensions(prhs[0]);
    ndims = mxGetNumberOfDimensions(prhs[0]);
    
    const mwSize *size_I = mxGetDimensions(prhs[0]);
    /////////////////////////////////////////////////
    mexPrintf("n elements = %d.\n", num_elements);
    
    /////////////////////////////////////////////////
    plhs[0] = mxCreateNumericArray(ndims, input_size, mxDOUBLE_CLASS, mxREAL);
    L = mxGetPr(plhs[0]);
    
    
    switch (class_id) {
        case mxLOGICAL_CLASS:
            compute_watershed((mxLogical *)I, (double *) mxGetData(prhs[1]),  num_elements, L, size_I);
            break;
            
        case mxUINT8_CLASS:
            compute_watershed((uint8_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        case mxUINT16_CLASS:
            compute_watershed((uint16_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        case mxUINT32_CLASS:
            compute_watershed((uint32_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        case mxINT8_CLASS:
            compute_watershed((int8_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        case mxINT16_CLASS:
            compute_watershed((int16_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        case mxINT32_CLASS:
            compute_watershed((int32_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        case mxSINGLE_CLASS:
            do_nan_check((float *)I, num_elements);
            compute_watershed((float *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        case mxDOUBLE_CLASS:
            do_nan_check((double *)I, num_elements);
            compute_watershed((double *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I);
            break;
            
        default:
            mxAssert(false, "Unexpected mxClassID in switch statement");
            break;
    }
}
