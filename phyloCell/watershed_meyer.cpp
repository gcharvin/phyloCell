// Copyright 2006-2007 The MathWorks, Inc.

#include "mex.h"
#include "neighborhood.h"
#include "FifoPriorityQueue.h"

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
inline void do_nan_check(_T *F,mwSize num_elements)
{
    for (mwSize p = 0; p < num_elements; p++)
    {
        if (mxIsNaN(F[p]))
        {
            mexErrMsgIdAndTxt("Images:watershed:expectedNonNaN",
                              "%s",
                              "Input image may not contain NaNs.");
        }
    }
}     

///////////////////////////////////////////////////////////////////////////////
// Algorithm reference: F. Meyer, "Topographic distance and watershed lines,"
// Signal Processing, 38:113-125, 1994.  Implemented using the description
// in R. Beare and G. Lehmann, "The watershed transform in ITK - discussion
// and new developments," The Insight Journal, 2006 January - June,
// www.insight-journal.org, downloaded 7-Jul-2006.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
void compute_watershed(_T *I, double *M, mwSize N, NeighborhoodWalker_T walker,
                       double *L)
{
    FifoPriorityQueue<mwSize, _T> queue(FifoPriorityItemCompareFcn<mwSize, _T>::LowestPriorityFirst);

    bool *S = new bool[N];
    for (mwSize p = 0; p < N; p++)
    {
        S[p] = false;
    }
    
    const double WSHED = 0.0;

    for (mwSize p = 0; p < N; p++)
    {
        L[p] = M[p];
        if (M[p] != WSHED)
        {
            mwSize q;
            S[p] = true;
            nhSetWalkerLocation(walker, p);
            while (nhGetNextInboundsNeighbor(walker, &q, NULL))
            {
                if ( (! S[q]) && (M[q] == WSHED) )
                {
                    S[q] = true;
                    queue.push(q, I[q]);
                }
            }
        }
    }

    while (! queue.isEmpty() )
    {
        mwSize q;
        mwSize p = queue.topData();
        _T  v = queue.topPriority();
        queue.pop();

        double label = WSHED;
        bool watershed = false;

        nhSetWalkerLocation(walker, p);
        while (nhGetNextInboundsNeighbor(walker, &q, NULL))
        {
            if ((L[q] != WSHED) && !watershed)
            {
                if ((label != WSHED) && (L[q] != label))
                {
                    watershed = true;
                }
                else
                {
                    label = L[q];
                }
            }
        }

        if (!watershed)
        {
            L[p] = label;
            nhSetWalkerLocation(walker, p);
            while (nhGetNextInboundsNeighbor(walker, &q, NULL))
            {
                if (!S[q])
                {
                    S[q] = true;
                    queue.push(q, std::max(I[q], v));
                }
            }
        }
    }

    delete[] S;
}

#define SIZE_MISMATCH_ID  "Images:watershed_meyer:sizeMismatch"
#define SIZE_MISMATCH_MSG "I must be the same size as L."

//////////////////////////////////////////////////////////////////////////////
// Do some input validation.  Validation of the type of I, as well as
// the connectivity argument, is done in toolbox/images/images/watershed.m.
//////////////////////////////////////////////////////////////////////////////
void check_inputs(int nrhs, const mxArray *prhs[])
{
    mwSize num_dims;

    if (nrhs != 3)
    {
        mexErrMsgIdAndTxt("Images:watershed_meyer:invalidNumInputs",
                          "%s",
                          "WATERSHED_MEYER needs 3 input arguments.");
    }

    if (!mxIsDouble(prhs[2]))
    {
        mexErrMsgIdAndTxt("Images:watershed_meyer:nondoubleMinimaLabels",
                          "%s",
                          "Third input argument to WATERSHED_MEYER must be double.");
    }

    num_dims = mxGetNumberOfDimensions(prhs[0]);
    if (num_dims != mxGetNumberOfDimensions(prhs[2]))
    {
        mexErrMsgIdAndTxt(SIZE_MISMATCH_ID, "%s", SIZE_MISMATCH_MSG);
    }

    const mwSize *size_I = mxGetDimensions(prhs[0]);
    const mwSize *size_L = mxGetDimensions(prhs[2]);
    for (mwSize k = 0; k < num_dims; k++)
    {
        if (size_I[k] != size_L[k])
        {
            mexErrMsgIdAndTxt(SIZE_MISMATCH_ID, "%s", SIZE_MISMATCH_MSG);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    void *I;
    double *L;
    mwSize num_elements;
    const mwSize *input_size;
    mwSize ndims;
    mxClassID class_id;
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;

    (void) nlhs;  /* unused parameter */

    check_inputs(nrhs, prhs);

    if (mxIsLogical(prhs[0]))
    {
        I = mxGetLogicals(prhs[0]);
    }
    else
    {
        I = mxGetData(prhs[0]);
    }
    num_elements = mxGetNumberOfElements(prhs[0]);
    class_id = mxGetClassID(prhs[0]);
    input_size = mxGetDimensions(prhs[0]);
    ndims = mxGetNumberOfDimensions(prhs[0]);
    
    plhs[0] = mxCreateNumericArray(ndims, input_size, mxDOUBLE_CLASS, mxREAL);
    L = mxGetPr(plhs[0]);

    nhood = nhMakeNeighborhood(prhs[1],NH_CENTER_MIDDLE_ROUNDDOWN);
    walker = nhMakeNeighborhoodWalker(nhood,input_size,ndims,NH_SKIP_CENTER);
    
    switch (class_id)
    {
    case mxLOGICAL_CLASS:
        compute_watershed((mxLogical *)I, (double *) mxGetData(prhs[2]),  num_elements,
                          walker, L);
        break;
        
    case mxUINT8_CLASS:
        compute_watershed((uint8_T *)I, (double *) mxGetData(prhs[2]), num_elements, walker, 
                          L);
        break;
        
    case mxUINT16_CLASS:
        compute_watershed((uint16_T *)I, (double *) mxGetData(prhs[2]), num_elements, 
                          walker, L);
        break;

    case mxUINT32_CLASS:
        compute_watershed((uint32_T *)I, (double *) mxGetData(prhs[2]), num_elements,
                          walker, L);
        break;

    case mxINT8_CLASS:
        compute_watershed((int8_T *)I, (double *) mxGetData(prhs[2]), num_elements, 
                          walker, L);
        break;

    case mxINT16_CLASS:
        compute_watershed((int16_T *)I, (double *) mxGetData(prhs[2]), num_elements,
                          walker, L);
        break;

    case mxINT32_CLASS:
        compute_watershed((int32_T *)I, (double *) mxGetData(prhs[2]), num_elements,
                          walker, L);
        break;

    case mxSINGLE_CLASS:
        do_nan_check((float *)I, num_elements);
        compute_watershed((float *)I, (double *) mxGetData(prhs[2]), num_elements, 
                          walker, L);
        break;

    case mxDOUBLE_CLASS:
        do_nan_check((double *)I, num_elements);
        compute_watershed((double *)I, (double *) mxGetData(prhs[2]), num_elements,
                          walker, L);
        break;

    default:
        mxAssert(false, "Unexpected mxClassID in switch statement");
        break;
    }

    nhDestroyNeighborhood(nhood);
    nhDestroyNeighborhoodWalker(walker);
}
