//#include "filter.h"


//void arm_biquad_cascade_df1_init_f32(
//  arm_biquad_casd_df1_inst_f32 * S,
//  unsigned char numStages,
//  float * pCoeffs,
//  float * pState)
//{
//  /* Assign filter stages */
//  S->numStages = numStages;

//  /* Assign coefficient pointer */
//  S->pCoeffs = pCoeffs;

//  /* Clear state buffer and size is always 4 * numStages */
//  memset(pState, 0, (4u * (unsigned) numStages) * sizeof(float));

//  /* Assign state pointer */
//  S->pState = pState;
//}

//void arm_biquad_cascade_df1_f32(
//  const arm_biquad_casd_df1_inst_f32 * S,
//  float * pSrc,
//  float * pDst,
//  unsigned blockSize)
//{
//  float *pIn = pSrc;                         /*  source pointer            */
//  float *pOut = pDst;                        /*  destination pointer       */
//  float *pState = S->pState;                 /*  pState pointer            */
//  float *pCoeffs = S->pCoeffs;               /*  coefficient pointer       */
//  float acc;                                 /*  Simulates the accumulator */
//  float b0, b1, b2, a1, a2;                  /*  Filter coefficients       */
//  float Xn1, Xn2, Yn1, Yn2;                  /*  Filter pState variables   */
//  float Xn;                                  /*  temporary input           */
//  unsigned sample, stage = S->numStages;         /*  loop counters             */



//  /* Run the below code for Cortex-M4 and Cortex-M3 */

//  do
//  {
//    /* Reading the coefficients */
//    b0 = *pCoeffs++;
//    b1 = *pCoeffs++;
//    b2 = *pCoeffs++;
//    a1 = *pCoeffs++;
//    a2 = *pCoeffs++;

//    /* Reading the pState values */
//    Xn1 = pState[0];
//    Xn2 = pState[1];
//    Yn1 = pState[2];
//    Yn2 = pState[3];

//    /* Apply loop unrolling and compute 4 output values simultaneously. */
//    /*      The variable acc hold output values that are being computed:
//     *
//     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
//     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
//     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
//     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
//     */

//    sample = blockSize >> 2u;

//    /* First part of the processing with loop unrolling.  Compute 4 outputs at a time.
//     ** a second loop below computes the remaining 1 to 3 samples. */
//    while(sample > 0u)
//    {
//      /* Read the first input */
//      Xn = *pIn++;

//      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
//      Yn2 = (b0 * Xn) + (b1 * Xn1) + (b2 * Xn2) + (a1 * Yn1) + (a2 * Yn2);

//      /* Store the result in the accumulator in the destination buffer. */
//      *pOut++ = Yn2;

//      /* Every time after the output is computed state should be updated. */
//      /* The states should be updated as:  */
//      /* Xn2 = Xn1    */
//      /* Xn1 = Xn     */
//      /* Yn2 = Yn1    */
//      /* Yn1 = acc   */

//      /* Read the second input */
//      Xn2 = *pIn++;

//      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
//      Yn1 = (b0 * Xn2) + (b1 * Xn) + (b2 * Xn1) + (a1 * Yn2) + (a2 * Yn1);

//      /* Store the result in the accumulator in the destination buffer. */
//      *pOut++ = Yn1;

//      /* Every time after the output is computed state should be updated. */
//      /* The states should be updated as:  */
//      /* Xn2 = Xn1    */
//      /* Xn1 = Xn     */
//      /* Yn2 = Yn1    */
//      /* Yn1 = acc   */

//      /* Read the third input */
//      Xn1 = *pIn++;

//      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
//      Yn2 = (b0 * Xn1) + (b1 * Xn2) + (b2 * Xn) + (a1 * Yn1) + (a2 * Yn2);

//      /* Store the result in the accumulator in the destination buffer. */
//      *pOut++ = Yn2;

//      /* Every time after the output is computed state should be updated. */
//      /* The states should be updated as: */
//      /* Xn2 = Xn1    */
//      /* Xn1 = Xn     */
//      /* Yn2 = Yn1    */
//      /* Yn1 = acc   */

//      /* Read the forth input */
//      Xn = *pIn++;

//      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
//      Yn1 = (b0 * Xn) + (b1 * Xn1) + (b2 * Xn2) + (a1 * Yn2) + (a2 * Yn1);

//      /* Store the result in the accumulator in the destination buffer. */
//      *pOut++ = Yn1;

//      /* Every time after the output is computed state should be updated. */
//      /* The states should be updated as:  */
//      /* Xn2 = Xn1    */
//      /* Xn1 = Xn     */
//      /* Yn2 = Yn1    */
//      /* Yn1 = acc   */
//      Xn2 = Xn1;
//      Xn1 = Xn;

//      /* decrement the loop counter */
//      sample--;

//    }

//    /* If the blockSize is not a multiple of 4, compute any remaining output samples here.
//     ** No loop unrolling is used. */
//    sample = blockSize & 0x3u;

//    while(sample > 0u)
//    {
//      /* Read the input */
//      Xn = *pIn++;

//      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
//      acc = (b0 * Xn) + (b1 * Xn1) + (b2 * Xn2) + (a1 * Yn1) + (a2 * Yn2);

//      /* Store the result in the accumulator in the destination buffer. */
//      *pOut++ = acc;

//      /* Every time after the output is computed state should be updated. */
//      /* The states should be updated as:    */
//      /* Xn2 = Xn1    */
//      /* Xn1 = Xn     */
//      /* Yn2 = Yn1    */
//      /* Yn1 = acc   */
//      Xn2 = Xn1;
//      Xn1 = Xn;
//      Yn2 = Yn1;
//      Yn1 = acc;

//      /* decrement the loop counter */
//      sample--;

//    }

//    /*  Store the updated state variables back into the pState array */
//    *pState++ = Xn1;
//    *pState++ = Xn2;
//    *pState++ = Yn1;
//    *pState++ = Yn2;

//    /*  The first stage goes from the input buffer to the output buffer. */
//    /*  Subsequent numStages  occur in-place in the output buffer */
//    pIn = pDst;

//    /* Reset the output pointer */
//    pOut = pDst;

//    /* decrement the loop counter */
//    stage--;

//  } while(stage > 0u);
//}
