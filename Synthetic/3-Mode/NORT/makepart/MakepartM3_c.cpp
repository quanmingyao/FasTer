#include "mex.h"

// #include "blas.h"
// #include <cstddef>
// using namespace::std;

void mexFunction(int nargout, mxArray *pargout[],
	int nargin, const mxArray *pargin[])
	// (subs, U, V, szTen, vals(output))
{
	// tensor subscript
	const double* pSubs = mxGetPr(pargin[0]);

	// input matrix
	const double* pU = mxGetPr(pargin[1]);
	const double* pV = mxGetPr(pargin[2]);

	const size_t iK = mxGetN(pargin[1]);
	// mexPrintf("iK:%d \n", iK);

	// tensor size
	const double* pTenSz = mxGetPr(pargin[3]);

	const size_t iI1 = (size_t)pTenSz[0];
	const size_t iI2 = (size_t)pTenSz[1];
	const size_t iI3 = (size_t)pTenSz[2];

	// partial tensor (output)
	double* pVals = mxGetPr(pargin[4]);

	// subscript size
	const size_t iMode = mxGetN(pargin[0]);
	const size_t iNnz = mxGetM(pargin[0]);

	// main loop
	const double* pISub1 = pSubs;
	const double* pISub2 = pISub1 + iNnz;
	const double* pISub3 = pISub2 + iNnz;
	double* ipVals = pVals;

	for (size_t i = 0; i != iNnz; i++)
	{
		// subscript for each mode
		const size_t iSub1 = (size_t)(*pISub1);
		const size_t iSub2 = (size_t)(*pISub2);
		const size_t iSub3 = (size_t)(*pISub3);

		const double* pui = pU + (iSub3 - 1);
		const double* pvi = pV + (iSub2 - 1)*iI1 + (iSub1 - 1);

		double iVai = 0;
		for (size_t k = 0; k != iK ; k++)
		{
			iVai = iVai + pui[k*iI3] * pvi[k*(iI1*iI2)];
		}
		*ipVals = iVai;

		// mexPrintf("%d, %d, %d \n", iSub1, iSub2, iSub3);

		pISub1++;
		pISub2++;
		pISub3++;
		ipVals++;
	}
}