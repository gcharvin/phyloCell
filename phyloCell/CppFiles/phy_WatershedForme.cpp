// Copyright 2006-2007 The MathWorks, Inc.
#include "mex.h"
#include "FifoPriorityQueue.h"
#include "matrix.h"
#include <math.h>

#define getMax(x, y)     ((x)>(y)?(x):(y))
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
        void compute_watershed(_T *I, double *M, mwSize N, double *L, const mwSize *size_I, double *outputTest, const double *numberPixels) {
    
    FifoPriorityQueue<mwSize, _T> queue(FifoPriorityItemCompareFcn<mwSize, _T>::LowestPriorityFirst);
    /*=================================================================
     *declaration initialisation */
    int np=int(numberPixels[0]); /*number of pixels we aproximate the elipse*/
    _T *I2; /*copie of the initial image*/
    I2=new _T[N];
    
    mwSize s1= size_I[0];/*first size of the image (y)*/
    mwSize s2= size_I[1];/*second size of the image (x)*/
    mwSize neigh[8];/*array of neighbour pixels*/
    mwSize edge;/*pixels on the edge will not be analysed*/
    mwSize p;/*index of actual pixel in image*/
    mwSize count=0;
    bool haveNeigh=false;
    mwSize *penalized; /*the penalized pixels*/
    mwSize *deleted; /*the pixels reanalized (deleted)*/
    mwSize *F; /*already treted pixels (=L without the pixels in queue)*/
    mwSize *countPixelsObj; /*the nombre of pixels in each region*/
    mwSize numObj=0; /*the number of regions*/
    mwSize obj=0;/*actual object*/
    
    
    /*ellipse variable declaration*/
    mwSize *x;
    mwSize *y;/*x and y of the pixels of actual object*/
    mwSize *ind; /*linear index of all the pixels of actual object*/
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
    
    
    /*penalisation variables*/
    double maxI=-1000000.0; /*the maximum value of the height in current object*/
    double l=0;/*added lenght to the ellipse foci*/
    mwSize ne=600;/*number of points to construct an ellipse*/
    double xe=0.0;
    double ye=0.0;
    mwSize *xea;
    mwSize *yea;/*x and y of the constructed ellipse*/
    mwSize *indea;/*liniar indice of ellipe*/
    
    /*ellipse contour variables*/
    xea=new mwSize[ne];
    yea=new mwSize[ne];
    indea=new mwSize[ne];
    
    /*initialize with 0 all variables*/
    penalized=new mwSize[N];
    deleted=new mwSize[N];
    
    F=new mwSize[N];
    x=new mwSize[N];
    y=new mwSize[N];
    ind=new mwSize[N];
    
    for (mwSize p = 0; p < N; p++) {
        penalized[p] = 0;
        deleted[p]=0;
        F[p]=0;
        x[p]=0;
        y[p]=0;
        ind[p]=0;
        I2[p]=I[p];/*make a copy of the input image*/
    }
    
    /*put in queue first marqueurs points*/
    for (mwSize p = 0; p < N; p++) {
        L[p] = M[p];
        numObj=getMax(numObj, mwSize(L[p]));/*the maximum value in L is the number of objects numObj*/
        if (M[p] != 0) {
            queue.push(p, I2[p]); /*add marqueurs to queue*/
            
        }
    }
    /*initialize countPixelsObj with 0*/
    countPixelsObj=new mwSize[numObj];
    for (mwSize p = 0; p < numObj; p++) {
        countPixelsObj[p]=0;
    }
    /*====================================================================
     *main WHILE */
    while (! queue.isEmpty() ) {
        /*pop the new element*/
        p = queue.topData();
        _T  v = queue.topPriority();
        queue.pop();
        
        
        if (F[p]!=0){
            continue;
        }
        
        /*L[p] can be zero after we deleted pixel p*/
        if (L[p]==0){
            continue;
        }
        
        /*re put in queue the penalized pixels, penalize xx times the same pixel*/
        if ((penalized[p]>=1) && (penalized[p]<=10)){
            if ((penalized[p]%2)==1){
                queue.push(p, I2[p]);
                penalized[p]++;
                continue;
            }
        }
        
        /*==============================================================
         * push the new neighbour pixels in the queue*/
        haveNeigh=false;
        edge=p%s1;
        if ((p>=s1) && (p<=N-s1) && ((edge!=0) && (edge!=s1-1))) { /*edge checking*/
            
            neigh[0]=p-s1;/*left*/
            neigh[1]=p-1;/*up*/
            neigh[2]=p+1;/*down*/
            neigh[3]=p+s1;/*right*/
            /*uncoment folowings for eight neigh*/
            //neigh[4]=p-s1-1;
            //neigh[5]=p-s1+1;
            //neigh[6]=p+s1-1;
            //neigh[7]=p+s1+1;
            haveNeigh=true;
        }
        
        /*if have neighbour, push them in que*/
        if (haveNeigh){
            for (mwSize j = 0; j < 4; j++) {  /*j<8 for eight neigh*/
                if (L[neigh[j]]==0) {
                    L[neigh[j]]=L[p];
                    queue.push(neigh[j], I2[neigh[j]]);
                }
            }
        }/*end if push neighbourg
         * ---------------------------------------------------------------*/
        
        
        F[p]=mwSize(L[p]);/*already visited*/
        
        obj=mwSize(L[p]);/*obj is the actual object*/
        
        countPixelsObj[obj-1]++;/*incremet count of pixels of each object */
        
        
        if ((L[p]>1)&&(countPixelsObj[obj-1]>=np)&&(countPixelsObj[obj-1]%np==0)){
            /*if not background(>1) and at each the numberPixels added to a region*/
            /*============================================================
             *calcul of ellipse parametres*/
            
            count=0;
            xbar=0;
            ybar=0;
            maxI=-1000000;
            
            /*analyse all pixels of image and, for each pixel having the same value as obj, kepp it in an array*/
            for (mwSize i = 0; i < N; i++){
                if (L[i]==obj){
                    x[count]=i/s1;
                    y[count]=i%s1;/*transform from linear index to x/ y index*/
                    ind[count]=i;
                    xbar+=x[count];
                    ybar+=y[count];
                    maxI=getMax(maxI, I2[p]);/*for the penalisation*/
                    count++;
                }
            } /*get x and y of the pixels from the actual object and their mean*/
            
            xbar/=double(count);
            ybar/=double(count);
            stats.Centroid[0]=xbar;
            stats.Centroid[1]=ybar;/*mean x,y*/
            /* Calculate normalized second central moments for the region. 1/12 is
             * the normalized second central moment of a pixel with unit length.*/
            uxx=0;
            uyy=0;
            uxy=0;
            for (mwSize i = 0; i < count; i++){
                difx = x[i] - xbar;
                dify = -(y[i] - ybar);
                /*This is negative for the orientation calculation(measured in the counter-clockwise direction).*/
                uxx += difx*difx;
                uyy += dify*dify;
                uxy += difx*dify;
            }
            
            uxx/=count;
            uxx+=1.0/12;
            uyy/=count;
            uyy+=1.0/12;
            uxy/=count;
            
            /* Calculate major axis length, minor axis length, and eccentricity.*/
            common = sqrt(pow((uxx - uyy), 2) + 4*uxy*uxy);
            stats.MajorAxisLength = 2*sqrt(2.0)*sqrt(uxx + uyy + common);
            stats.MinorAxisLength = 2*sqrt(2.0)*sqrt(uxx + uyy - common);
            
            stats.Eccentricity = 2*sqrt(pow((stats.MajorAxisLength/2), 2)-pow((stats.MinorAxisLength/2), 2))/stats.MajorAxisLength;
            
            /* Calculate orientation.*/
            if (uyy > uxx){
                num = uyy - uxx + common;
                den = 2*uxy;
            }
            else{
                num = 2*uxy;
                den = uxx - uyy + common;
            }
            
            if ((num == 0) && (den == 0))
                stats.Orientation = 0;
            else
                stats.Orientation = atan(num/den);
            
            
            /*calculate foci of ellipse f1 and f2 (f1[O] is f1x, f1[1] is f1y)*/
            A=stats.MajorAxisLength*1;
            a=stats.MinorAxisLength*1;
            alpha=stats.Orientation;
            
            stats.f1[0]=-sqrt(A*A-a*a)/2*cos(alpha)+stats.Centroid[0];
            stats.f1[1]=sqrt(A*A-a*a)/2*sin(alpha)+stats.Centroid[1];
            stats.f2[0]=sqrt(A*A-a*a)/2*cos(alpha)+stats.Centroid[0];
            stats.f2[1]=-sqrt(A*A-a*a)/2*sin(alpha)+stats.Centroid[1];
            
            /*============================================================
             *end of calcul of ellipse parametres*/
            
            
            /*calcul de coefficients de penalisation
             *=====================================================*/
            double cosAlpha;
            double sinAlpha;
            bool colision=false;
            cosAlpha=cos(alpha);
            sinAlpha=sin(alpha);
            /*cntruct the ellipse contour with ne pouints (600)*/
            for (mwSize i = 0; i < ne; i++){
                //x, y ellipse
                xe=A/2*cos((2*pi*i/ne));
                ye=a/2*sin((2*pi*i/ne));
                
                //x,y ellipse after rotation
                xea[i]=floor(xe*cosAlpha+ye*sinAlpha+stats.Centroid[0]+0.5);/*(round to the clossest integer)*/
                yea[i]=floor(-xe*sinAlpha+ye*cosAlpha+stats.Centroid[1]+0.5);
                
                /*check the validity of x and y (inside image)*/
                if (xea[i]<0)
                    xea[i]=0;
                if (yea[i]<0)
                    yea[i]=0;
                if (xea[i]>=s2)
                    xea[i]=s2-1;
                if (yea[i]>=s1)
                    yea[i]=s1-1;
                
                indea[i]=xea[i]*s1+yea[i]; //calculate the linear index
                maxI=getMax(maxI, I2[indea[i]]); //panalisation is the maximum on the ellipse contour
                
                /*test if ellipse collision*/
                if ((L[indea[i]]!=0)&&(L[indea[i]]!=obj))
                    colision=true;
            }
            
            /*begin to penalize only if the ellipse aproximation intertsect with others objects*/
            if (colision){
                /*===========================================================
                 *Penalize the pixels*/
                
                /*check all tha pixels of the actual object*/
                for (mwSize i = 0; i < count; i++){
                    /*calculate the summ of distance to the elipse foci*/
                    l=sqrt(pow((x[i]-stats.f1[0]), 2)+pow((y[i]-stats.f1[1]), 2))+sqrt(pow((x[i]-stats.f2[0]), 2)+pow((y[i]-stats.f2[1]), 2));
                    p=ind[i];
                    
                    /*if on the contour*/
                    if (l>=(A-1)&& l<=(A+1)){
                        penalized[p]++;
                        
                        //if (penalized[p]>2)
                        //maxI=maxI-0.1;
                        
                        if (penalized[p]<=10)
                            I2[p]=maxI; /*penalize the pixel*/
                        
                        
                        if ((F[p]==obj) && (penalized[p]<=10)){
                            
                            I2[p]=maxI; /*penalize the pixel*/
                            queue.push(p, I2[p]);
                            F[p]=0;
                            penalized[p]++;
                        }
                    }

                    /*if outside the contour*/
                    if (l>(A+1)){
                        if (deleted[p]<=10){
                            deleted[p]++;
                            
                            /*for each deleted pixel check if the neighbours were already treated*/
                            haveNeigh=false;
                            edge=p%s1;
                            if ((p>=s1) && (p<=N-s1) && ((edge!=0) && (edge!=s1-1))) {
                                
                                neigh[0]=p-s1;/*left*/
                                neigh[1]=p-1;/*up*/
                                neigh[2]=p+1;/*down*/
                                neigh[3]=p+s1;/*right*/
                                haveNeigh=true;
                            }
                            
                            if (haveNeigh){
                                for (mwSize j = 0; j < 4; j++) {  /*j<8 for eight neigh*/
                                    if ((F[neigh[j]]!=0)&&(F[neigh[j]]!=obj)) {/*if already treated an not in actual object*/
                                        F[neigh[j]]=0;
                                        outputTest[neigh[j]]++;
                                        queue.push(neigh[j], I2[neigh[j]]);/*retreat those pixes*/
                                    }
                                }
                            }
                            
                            //I2[p]=maxI; /*penalize the pixel*/
                            L[p]=0;
                            F[p]=0;
                            //countPixelsObj[obj-1]--;
                        }
                    }
                    if (l<(A-1)&&(F[p]==0) && (penalized[p]>=1)){
                        /*if inside the contour*/
                        queue.push(p, I[p]); /*if penalized pixels depassed then their priority is the one non penalized*/
                    }
                    
                }
                /*===========================================================
                 *end penalize pixels*/
            }/*end if colision*/
            
        } /*end if approximarion*/
        
        
        
        
        
    }/*end main while*/
    
    
    for (mwSize i = 0; i < N; i++){
        //outputTest[i]=deleted[i];
        // mexPrintf("%d,%d,%d;\n",xea[i],yea[i],indea[i]);
    }
    //mexPrintf("test%d \n",mwSize(floor(5.5673+0.5)));
    
    //mxFree(x);
    //mxFree(y);
    delete[] y;
    delete[] x;
    delete[] penalized;
    delete[] deleted;
    delete[] F;
    delete[] countPixelsObj;
    delete[] ind;
    delete[] xea;
    delete[] yea;
    delete[] indea;
    delete[] I2;
    
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
    double *outTest;
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
    
    const double *numberPixels=mxGetPr(prhs[2]);
    /////////////////////////////////////////////////
    mexPrintf("n elements = %d.\n", num_elements);
    
    
    /////////////////////////////////////////////////
    plhs[0] = mxCreateNumericArray(ndims, input_size, mxDOUBLE_CLASS, mxREAL);
    L = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(ndims, input_size, mxDOUBLE_CLASS, mxREAL);
    outTest = mxGetPr(plhs[1]);
    
    switch (class_id) {
        case mxLOGICAL_CLASS:
            compute_watershed((mxLogical *)I, (double *) mxGetData(prhs[1]),  num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxUINT8_CLASS:
            compute_watershed((uint8_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxUINT16_CLASS:
            compute_watershed((uint16_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxUINT32_CLASS:
            compute_watershed((uint32_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxINT8_CLASS:
            compute_watershed((int8_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxINT16_CLASS:
            compute_watershed((int16_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxINT32_CLASS:
            compute_watershed((int32_T *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxSINGLE_CLASS:
            do_nan_check((float *)I, num_elements);
            compute_watershed((float *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        case mxDOUBLE_CLASS:
            do_nan_check((double *)I, num_elements);
            compute_watershed((double *)I, (double *) mxGetData(prhs[1]), num_elements, L, size_I, outTest, numberPixels);
            break;
            
        default:
            mxAssert(false, "Unexpected mxClassID in switch statement");
            break;
    }
}
