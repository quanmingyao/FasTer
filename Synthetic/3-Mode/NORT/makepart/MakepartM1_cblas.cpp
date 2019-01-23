#include "mex.h"

#include "blas.h"
#include <cstddef>
using namespace::std;

void mexFunction(int nargout, mxArray *pargout[],
	int nargin, const mxArray *pargin[])
	// (subs, U, V, szTen, vals(output))
{
	// tensor subscript
	const double* pSubs = mxGetPr(pargin[0]);

	// input matrix
	const double* pU = mxGetPr(pargin[1]);
	const double* pV = mxGetPr(pargin[2]);

	const ptrdiff_t iK = mxGetN(pargin[1]);
	// mexPrintf("iK:%d \n", iK);

	// tensor size
	const double* pTenSz = mxGetPr(pargin[3]);

	const ptrdiff_t iI1 = (ptrdiff_t)pTenSz[0];
	const ptrdiff_t iI2 = (ptrdiff_t)pTenSz[1];
	const ptrdiff_t iI3 = (ptrdiff_t)pTenSz[2];

	// partial tensor (output)
	double* pVals = mxGetPr(pargin[4]);

	// subscript size
	const ptrdiff_t iMode = mxGetN(pargin[0]);
	const ptrdiff_t iNnz = mxGetM(pargin[0]);

	// main loop
	const double* pISub1 = pSubs;
	const double* pISub2 = pISub1 + iNnz;
	const double* pISub3 = pISub2 + iNnz;
	double* ipVals = pVals;

	for (ptrdiff_t i = 0; i != iNnz; i++)
	{
		// subscript for each mode
		const ptrdiff_t iSub1 = (ptrdiff_t)(*pISub1);
		const ptrdiff_t iSub2 = (ptrdiff_t)(*pISub2);
		const ptrdiff_t iSub3 = (ptrdiff_t)(*pISub3);

		const double* pui = pU + (iSub1 - 1);
		const double* pvi = pV + (iSub3 - 1)*iI2 + (iSub2 - 1);

		/*double iVai = 0;
		for (ptrdiff_t k = 0; k != iK ; k++)
		{
			iVai = iVai + pui[k*iI1] * pvi[k*(iI2*iI3)];
		}
		*ipVals = iVai;*/

		const ptrdiff_t incy = iI2*iI3;
		*ipVals = ddot(&iK, pui, &iI1, pvi, &incy);

		// mexPrintf("%d, %d, %d \n", iSub1, iSub2, iSub3);

		pISub1++;
		pISub2++;
		pISub3++;
		ipVals++;
	}
}