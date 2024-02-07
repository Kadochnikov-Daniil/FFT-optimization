#include "dft.h"

extern void dft8Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    //dft4 x dft2
    //DFT2 * 4
    cfloat32_t firstTmpDst[8];
    firstTmpDst[0].re = pSrc[0].re + pSrc[4].re;
    firstTmpDst[0].im = pSrc[0].im + pSrc[4].im;
    firstTmpDst[1].re = pSrc[0].re - pSrc[4].re;
    firstTmpDst[1].im = pSrc[0].im - pSrc[4].im;

    firstTmpDst[2].re = pSrc[1].re + pSrc[5].re;
    firstTmpDst[2].im = pSrc[1].im + pSrc[5].im;
    firstTmpDst[3].re = pSrc[1].re - pSrc[5].re;
    firstTmpDst[3].im = pSrc[1].im - pSrc[5].im;

    firstTmpDst[4].re = pSrc[2].re + pSrc[6].re;
    firstTmpDst[4].im = pSrc[2].im + pSrc[6].im;
    firstTmpDst[5].re = pSrc[2].im - pSrc[6].im;
    firstTmpDst[5].im = pSrc[6].re - pSrc[2].re;

    firstTmpDst[6].re = pSrc[3].re + pSrc[7].re;
    firstTmpDst[6].im = pSrc[3].im + pSrc[7].im;
    firstTmpDst[7].re = pSrc[7].re - pSrc[3].re;
    firstTmpDst[7].im = pSrc[7].im - pSrc[3].im;

    //coefficient multiplication
    float pi_4 = 0.707106781186548;
    float im = (firstTmpDst[3].im - firstTmpDst[3].re) * pi_4;
    firstTmpDst[3].re = (firstTmpDst[3].re + firstTmpDst[3].im) * pi_4;
    firstTmpDst[3].im = im;

    im = (firstTmpDst[7].im + firstTmpDst[7].re) * pi_4;
    firstTmpDst[7].re = (firstTmpDst[7].re - firstTmpDst[7].im) * pi_4;
    firstTmpDst[7].im = im;

    //DFT4 * 2
    cfloat32_t secondTmpDst[8];
    secondTmpDst[0].re = firstTmpDst[0].re + firstTmpDst[4].re;
    secondTmpDst[0].im = firstTmpDst[0].im + firstTmpDst[4].im;
    secondTmpDst[1].re = firstTmpDst[0].re - firstTmpDst[4].re;
    secondTmpDst[1].im = firstTmpDst[0].im - firstTmpDst[4].im;
    secondTmpDst[2].re = firstTmpDst[2].re + firstTmpDst[6].re;
    secondTmpDst[2].im = firstTmpDst[2].im + firstTmpDst[6].im;
    secondTmpDst[3].re = firstTmpDst[2].re - firstTmpDst[6].re;
    secondTmpDst[3].im = firstTmpDst[2].im - firstTmpDst[6].im;

    secondTmpDst[4].re = firstTmpDst[1].re + firstTmpDst[5].re;
    secondTmpDst[4].im = firstTmpDst[1].im + firstTmpDst[5].im;
    secondTmpDst[5].re = firstTmpDst[1].re - firstTmpDst[5].re;
    secondTmpDst[5].im = firstTmpDst[1].im - firstTmpDst[5].im;
    secondTmpDst[6].re = firstTmpDst[3].re + firstTmpDst[7].re;
    secondTmpDst[6].im = firstTmpDst[3].im + firstTmpDst[7].im;
    secondTmpDst[7].re = firstTmpDst[3].re - firstTmpDst[7].re;
    secondTmpDst[7].im = firstTmpDst[3].im - firstTmpDst[7].im;

    pDst[0].re = secondTmpDst[0].re + secondTmpDst[2].re;
    pDst[0].im = secondTmpDst[0].im + secondTmpDst[2].im;
    pDst[2].re = secondTmpDst[1].re + secondTmpDst[3].im;
    pDst[2].im = secondTmpDst[1].im - secondTmpDst[3].re;
    pDst[4].re = secondTmpDst[0].re - secondTmpDst[2].re;
    pDst[4].im = secondTmpDst[0].im - secondTmpDst[2].im;
    pDst[6].re = secondTmpDst[1].re - secondTmpDst[3].im;
    pDst[6].im = secondTmpDst[1].im + secondTmpDst[3].re;

    pDst[1].re = secondTmpDst[4].re + secondTmpDst[6].re;
    pDst[1].im = secondTmpDst[4].im + secondTmpDst[6].im;
    pDst[3].re = secondTmpDst[5].re + secondTmpDst[7].im;
    pDst[3].im = secondTmpDst[5].im - secondTmpDst[7].re;
    pDst[5].re = secondTmpDst[4].re - secondTmpDst[6].re;
    pDst[5].im = secondTmpDst[4].im - secondTmpDst[6].im;
    pDst[7].re = secondTmpDst[5].re - secondTmpDst[7].im;
    pDst[7].im = secondTmpDst[5].im + secondTmpDst[7].re;
    return;
}