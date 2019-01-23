#include "mex.h"
#include "blas.h"
#include <cstddef>

using namespace::std;

/* The gateway function */
void mexFunction(int nargout, mxArray *pargout[],
                 int nargin , const mxArray *pargin[])
// (subs, vals, U, szTen)
{
	// tensor subscript
    const double* pSubs = mxGetPr(pargin[0]);

	// partial tensor
    double* pVals = mxGetPr(pargin[1]);
    
    // input matrix
    double* pU = mxGetPr(pargin[2]);

	const ptrdiff_t iK = mxGetN(pargin[2]);

	// mexPrintf("iK:%d \n", iK);
    
    // tensor size
    double* pTenSz = mxGetPr(pargin[3]);

	const ptrdiff_t iI1 = (ptrdiff_t) pTenSz[0];
	const ptrdiff_t iI2 = (ptrdiff_t) pTenSz[1];
	const ptrdiff_t iI3 = (ptrdiff_t) pTenSz[2];
    
    // output matrix
	pargout[0] = mxCreateDoubleMatrix((mwSize)iK, (mwSize)(iI2*iI3), mxREAL);
    double* pX = mxGetPr(pargout[0] );
    
	// subscript size
    const size_t iMode = mxGetN(pargin[0]);
	const size_t iNnz  = mxGetM(pargin[0]);


	const double* pISub1 = pSubs;
	const double* pISub2 = pISub1 + iNnz;
	const double* pISub3 = pISub2 + iNnz;
	const double* ipVals = pVals;

	ptrdiff_t incx = 1;

	for (size_t i = 0; i != iNnz; i++)
	{
		// subscript for each mode
		const size_t iSub1 = (size_t) (*pISub1);
		const size_t iSub2 = (size_t) (*pISub2);
		const size_t iSub3 = (size_t) (*pISub3);

		// mexPrintf("%d, %d, %d \n", iSub1, iSub2, iSub3);

		const size_t ipi = (iSub2 - 1) + iI2*(iSub3 - 1);

		// const double dVali = *ipVals;

		double* pxi = pX + iK*ipi;
		double* pui = pU + (iSub1 - 1);

		/* for (size_t k = 0; k < iK; k++)
		{
			pxi[k] = pxi[k] + dVali * pui[k*iI1];

			// mexPrintf("%.4f ", pui[k*iI1]);
		}
		// mexPrintf("\n"); */

		
		daxpy(&iK, ipVals, pui, &iI1, pxi, &incx);
		
		pISub1++;
		pISub2++;
		pISub3++;
		ipVals++;
	}
}