#include "mex.h"
#include "math.h"
#include "stdio.h"

#ifndef M_PI
    #define M_PI           3.14159265358979323846
#endif

// #define DEBUG        
//---------------------------------------------------------------------------
float am_fsqrt(float x)  // 0,00079 on sqrt(220.54) , t = 220, tref = 580
//http://www.codemaestro.com/reviews/review00000105.html
{
    float x_in = x ;
    float xhalf = 0.5f*x;
    int i = *(int*)&x; // get bits for floating value
    i = 0x5f375a86- (i>>1); // gives initial guess y0
    x = *(float*)&i; // convert bits back to float
    x = x*(1.5f-xhalf*x*x); // Newton step, repeating increases accuracy
    return x*x_in;
}
//---------------------------------------------------------------------------

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // parameters: mu, w, Cov, G = scalar, assumed also (F = 1)
    // mu --> [data1; data2; data3; data4; ... ]
    // Cov ---> [data2;data2;data3;...]
     
 double *muValues, *cValues, *wValues, *gValue, tmp, *ptr_out ;
 
 double A, dm, m, trFA_in, trFA_out, detA ;
 int i, j, rowLen, colLen, k_i, idx_i, idx_j ;
 double eta, I, const_norm ;
 
 if ( nrhs != 4 )
     mexErrMsgTxt("To few input parameters! (Mu, w, Cov)") ;
 
//     	n_elms = mxGetNumberOfElements (prhs[2]);
 
 // debug info
 #ifdef DEBUG
 printf("Number of inputs: %d\n", nrhs) ;
 for ( int i = 0 ; i < nrhs-1 ; ++i)
    printf("P1 [col, row]: [%d, %d]\n", mxGetN(prhs[i]), mxGetM(prhs[i])) ; 
 fflush(stdout) ;
 #endif
 
 
 //Get matrix mu
 muValues = mxGetPr(prhs[0]);
 colLen = (int)mxGetN(prhs[0]);
 rowLen = (int)mxGetM(prhs[0]);
  
 if ( mxGetM(prhs[0]) != mxGetM(prhs[2]) )
      mexErrMsgTxt("rows of Mu and Cov not equal!") ; 
 if ( mxGetN(prhs[0]) != mxGetN(prhs[1]) )
      mexErrMsgTxt("columns of Mu and Cov not equal!") ;
 if ( mxGetN(prhs[0]) != mxGetN(prhs[2]) )
      mexErrMsgTxt("columns of Mu and w not equal!") ;
 
 wValues = mxGetPr(prhs[1]);
 cValues = mxGetPr(prhs[2]);
 gValue = mxGetPr(prhs[3]);
 
 #ifdef DEBUG
 for(int i = 0 ; i < colLen ; i++)
     printf("%f, ", wValues[i]) ;
 printf("Starting for loop...") ;
 #endif
 
 const_norm = pow(1.0/(2.0*M_PI), (((double)rowLen)/2.0)) ;
 I = 0.0 ;
 for( i = 0 ; i < colLen ; i++ ) {   
     idx_i = i*rowLen ;
//      printf("Vector mean %d: %f, %f \n", i, muValues[idx_i + 0], muValues[idx_i + 1]) ;
     
     for ( j = i ; j < colLen ; j++ ) {
         idx_j = j*rowLen ;
         if ( i == j ) {
             eta = 1.0 ;
         } else {
             eta = 2.0 ;
         }
         
         m = 0.0 ;
         trFA_in = 0.0 ;
         trFA_out = 0.0 ;
         detA = 1.0 ;
         for ( k_i = 0 ; k_i < rowLen ; ++k_i ) {
             A = 1.0 / (gValue[0] + cValues[idx_i + k_i] + cValues[idx_j + k_i]) ;
     
             dm = muValues[idx_i + k_i] - muValues[idx_j + k_i] ;
             m += dm*dm*A ;
             
             trFA_in += A*A ;
             trFA_out += A ;
             detA *= A ;
         }
         trFA_out *= trFA_out ;
          
       
//         I += wValues[i] * wValues[j] * eta * ( const_norm*am_fsqrt(detA)*exp(-0.5*m) ) * ( 2*trFA_in*(1.0-2.0*m) + (1.0-m)*(1.0-m)*trFA_out ) ;
         I += wValues[i] * wValues[j] * eta * ( const_norm*sqrt(detA)*exp(-0.5*m) ) * ( 2*trFA_in*(1.0-2.0*m) + (1.0-m)*(1.0-m)*trFA_out ) ;


//         if ( i == 0 &&  j == 1 ) {    
//             printf("trFA_in %f\n",trFA_in) ;
//             printf("trFA_out %f\n",trFA_out) ;
//             printf("m %f\n",m) ;
//             printf("Value: %f\n", tmp) ;
//         }
             
     }
 }
 #ifdef DEBUG
    printf("End of for loop...") ;
 #endif
         
 if ( nlhs > 0 ) {
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  
    ptr_out = (double *) mxGetData(plhs[0]) ;
    ptr_out[0] = (double)I ;        
 } else {
    printf("\n Result: %f", I ) ;    
 }
 
/* if ( nlhs > 0 ) {
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  
    outArray = mxGetPr(plhs[0]);
    outArray[0] = I ;
 }*/
}
