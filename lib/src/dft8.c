#include "dft.h"

extern void dft8Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    refDftFwd(pSrc, pDst, 8);
}
