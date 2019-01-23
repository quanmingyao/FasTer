#include "mex.h"

void mexFunction(int nargout, mxArray *pargout[],
                 int nargin , const mxArray *pargin[])
// (subs, vals, U, szTen)
{
	// tensor subscript
    const double* pSubs = mxGetPr(pargin[0]);

	// partial tensor
    const double* pVals = mxGetPr(pargin[1]);
    
    // input matrix
    const double* pU = mxGetPr(pargin[2]);

	const size_t iK = mxGetN(pargin[2]);

	// mexPrintf("iK:%d \n", iK);
    
    // tensor size
    const double* pTenSz = mxGetPr(pargin[3]);

	const size_t iI1 = (size_t) pTenSz[0];
	const size_t iI2 = (size_t) pTenSz[1];
	const size_t iI3 = (size_t) pTenSz[2];
    
    // output matrix
	pargout[0] = mxCreateDoubleMatrix((mwSize)iK, (mwSize)(iI1*iI2), mxREAL);
    double* pX = mxGetPr(pargout[0] );
    
	// subscript size
    const size_t iMode = mxGetN(pargin[0]);
	const size_t iNnz  = mxGetM(pargin[0]);


	const double* pISub1 = pSubs;
	const double* pISub2 = pISub1 + iNnz;
	const double* pISub3 = pISub2 + iNnz;
	const double* ipVals = pVals;

	// main loop
	for (size_t i = 0; i != iNnz; i++)
	{
		// subscript for each mode
		const size_t iSub1 = (size_t) (*pISub1);
		const size_t iSub2 = (size_t) (*pISub2);
		const size_t iSub3 = (size_t) (*pISub3);

		// mexPrintf("%d, %d, %d \n", iSub1, iSub2, iSub3);

		const size_t ipi = (iSub1 - 1) + iI1*(iSub2 - 1);
		const double* pui = pU + (iSub3 - 1);
		double* pxi = pX + iK*ipi;
		
		for (size_t k = 0; k != iK; k++)
		{
			pxi[k] = pxi[k] + (*ipVals) * pui[k*iI3];

			// mexPrintf("%.4f ", pui[k*iI1]);
		}
		// mexPrintf("\n");
		
		pISub1++;
		pISub2++;
		pISub3++;
		ipVals++;
	}
}