#include "dft.h"

extern void dft16Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    refDftFwd(pSrc, pDst, 16);
}
