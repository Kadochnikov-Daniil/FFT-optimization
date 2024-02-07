#include "dft.h"

extern void dft16Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    //dft4 x dft4
    //DFT4 * 4
    cfloat32_t firstTmpDst[16];
    firstTmpDst[0].re = pSrc[0].re + pSrc[8].re;
    firstTmpDst[0].im = pSrc[0].im + pSrc[8].im;
    firstTmpDst[1].re = pSrc[0].re - pSrc[8].re;
    firstTmpDst[1].im = pSrc[0].im - pSrc[8].im;
    firstTmpDst[2].re = pSrc[4].re + pSrc[12].re;
    firstTmpDst[2].im = pSrc[4].im + pSrc[12].im;
    firstTmpDst[3].re = pSrc[4].re - pSrc[12].re;
    firstTmpDst[3].im = pSrc[4].im - pSrc[12].im;

    firstTmpDst[4].re = pSrc[1].re + pSrc[9].re;
    firstTmpDst[4].im = pSrc[1].im + pSrc[9].im;
    firstTmpDst[5].re = pSrc[1].re - pSrc[9].re;
    firstTmpDst[5].im = pSrc[1].im - pSrc[9].im;
    firstTmpDst[6].re = pSrc[5].re + pSrc[13].re;
    firstTmpDst[6].im = pSrc[5].im + pSrc[13].im;
    firstTmpDst[7].re = pSrc[5].re - pSrc[13].re;
    firstTmpDst[7].im = pSrc[5].im - pSrc[13].im;

    firstTmpDst[8].re = pSrc[2].re + pSrc[10].re;
    firstTmpDst[8].im = pSrc[2].im + pSrc[10].im;
    firstTmpDst[9].re = pSrc[2].re - pSrc[10].re;
    firstTmpDst[9].im = pSrc[2].im - pSrc[10].im;
    firstTmpDst[10].re = pSrc[6].re + pSrc[14].re;
    firstTmpDst[10].im = pSrc[6].im + pSrc[14].im;
    firstTmpDst[11].re = pSrc[6].re - pSrc[14].re;
    firstTmpDst[11].im = pSrc[6].im - pSrc[14].im;

    firstTmpDst[12].re = pSrc[3].re + pSrc[11].re;
    firstTmpDst[12].im = pSrc[3].im + pSrc[11].im;
    firstTmpDst[13].re = pSrc[3].re - pSrc[11].re;
    firstTmpDst[13].im = pSrc[3].im - pSrc[11].im;
    firstTmpDst[14].re = pSrc[7].re + pSrc[15].re;
    firstTmpDst[14].im = pSrc[7].im + pSrc[15].im;
    firstTmpDst[15].re = pSrc[15].re - pSrc[7].re;
    firstTmpDst[15].im = pSrc[7].im - pSrc[15].im;

    cfloat32_t tmpDst[16];
    tmpDst[0].re = firstTmpDst[0].re + firstTmpDst[2].re;
    tmpDst[0].im = firstTmpDst[0].im + firstTmpDst[2].im;
    tmpDst[1].re = firstTmpDst[1].re + firstTmpDst[3].im;
    tmpDst[1].im = firstTmpDst[1].im - firstTmpDst[3].re;
    tmpDst[2].re = firstTmpDst[0].re - firstTmpDst[2].re;
    tmpDst[2].im = firstTmpDst[0].im - firstTmpDst[2].im;
    tmpDst[3].re = firstTmpDst[1].re - firstTmpDst[3].im;
    tmpDst[3].im = firstTmpDst[1].im + firstTmpDst[3].re;

    tmpDst[4].re = firstTmpDst[4].re + firstTmpDst[6].re;
    tmpDst[4].im = firstTmpDst[4].im + firstTmpDst[6].im;
    tmpDst[5].re = firstTmpDst[5].re + firstTmpDst[7].im;
    tmpDst[5].im = firstTmpDst[5].im - firstTmpDst[7].re;
    tmpDst[6].re = firstTmpDst[4].re - firstTmpDst[6].re;
    tmpDst[6].im = firstTmpDst[4].im - firstTmpDst[6].im;
    tmpDst[7].re = firstTmpDst[5].re - firstTmpDst[7].im;
    tmpDst[7].im = firstTmpDst[5].im + firstTmpDst[7].re;

    tmpDst[8].re = firstTmpDst[8].re + firstTmpDst[10].re;
    tmpDst[8].im = firstTmpDst[8].im + firstTmpDst[10].im;
    tmpDst[9].re = firstTmpDst[9].re + firstTmpDst[11].im;
    tmpDst[9].im = firstTmpDst[9].im - firstTmpDst[11].re;
    tmpDst[10].re = firstTmpDst[8].im - firstTmpDst[10].im;
    tmpDst[10].im = firstTmpDst[10].re - firstTmpDst[8].re;
    tmpDst[11].re = firstTmpDst[11].im - firstTmpDst[9].re;
    tmpDst[11].im = firstTmpDst[9].im + firstTmpDst[11].re;

    tmpDst[12].re = firstTmpDst[12].re + firstTmpDst[14].re;
    tmpDst[12].im = firstTmpDst[12].im + firstTmpDst[14].im;
    tmpDst[13].re = firstTmpDst[13].re + firstTmpDst[15].im;
    tmpDst[13].im = firstTmpDst[13].im + firstTmpDst[15].re;
    tmpDst[14].re = firstTmpDst[14].re - firstTmpDst[12].re;
    tmpDst[14].im = firstTmpDst[12].im - firstTmpDst[14].im;
    tmpDst[15].re = firstTmpDst[13].re - firstTmpDst[15].im;
    tmpDst[15].im = firstTmpDst[15].re - firstTmpDst[13].im;

    //coefficient multiplication
    float pi_8 = 0.923879532511287;
    float pi_4 = 0.707106781186548;
    float pi3_8 = 0.38268343236509;

    float im = tmpDst[5].im * pi_8 - tmpDst[5].re * pi3_8;
    tmpDst[5].re = tmpDst[5].re * pi_8 + tmpDst[5].im * pi3_8;
    tmpDst[5].im = im;

    im = (tmpDst[6].im - tmpDst[6].re) * pi_4;
    tmpDst[6].re = (tmpDst[6].re + tmpDst[6].im) * pi_4;
    tmpDst[6].im = im;

    im = tmpDst[7].im * pi3_8 - tmpDst[7].re * pi_8;
    tmpDst[7].re = tmpDst[7].re * pi3_8 + tmpDst[7].im * pi_8;
    tmpDst[7].im = im;

    im = (tmpDst[9].im - tmpDst[9].re) * pi_4;
    tmpDst[9].re = (tmpDst[9].re + tmpDst[9].im) * pi_4;
    tmpDst[9].im = im;

    im = (tmpDst[11].re - tmpDst[11].im) * pi_4;
    tmpDst[11].re = (tmpDst[11].im + tmpDst[11].re) * pi_4;
    tmpDst[11].im = im;

    im = tmpDst[13].im * pi3_8 - tmpDst[13].re * pi_8;
    tmpDst[13].re = tmpDst[13].re * pi3_8 + tmpDst[13].im * pi_8;
    tmpDst[13].im = im;

    im = (tmpDst[14].re - tmpDst[14].im ) * pi_4;
    tmpDst[14].re = (tmpDst[14].re + tmpDst[14].im) * pi_4;
    tmpDst[14].im = im;

    im = tmpDst[15].re * pi3_8 + tmpDst[15].im * pi_8;
    tmpDst[15].re = tmpDst[15].im * pi3_8 - tmpDst[15].re * pi_8;
    tmpDst[15].im = im;

    //DFT4 * 4
    cfloat32_t secondTmpDst[16];
    secondTmpDst[0].re = tmpDst[0].re + tmpDst[8].re;
    secondTmpDst[0].im = tmpDst[0].im + tmpDst[8].im;
    secondTmpDst[1].re = tmpDst[0].re - tmpDst[8].re;
    secondTmpDst[1].im = tmpDst[0].im - tmpDst[8].im;
    secondTmpDst[2].re = tmpDst[4].re + tmpDst[12].re;
    secondTmpDst[2].im = tmpDst[4].im + tmpDst[12].im;
    secondTmpDst[3].re = tmpDst[4].re - tmpDst[12].re;
    secondTmpDst[3].im = tmpDst[4].im - tmpDst[12].im;

    secondTmpDst[4].re = tmpDst[1].re + tmpDst[9].re;
    secondTmpDst[4].im = tmpDst[1].im + tmpDst[9].im;
    secondTmpDst[5].re = tmpDst[1].re - tmpDst[9].re;
    secondTmpDst[5].im = tmpDst[1].im - tmpDst[9].im;
    secondTmpDst[6].re = tmpDst[5].re + tmpDst[13].re;
    secondTmpDst[6].im = tmpDst[5].im + tmpDst[13].im;
    secondTmpDst[7].re = tmpDst[5].re - tmpDst[13].re;
    secondTmpDst[7].im = tmpDst[5].im - tmpDst[13].im;

    secondTmpDst[8].re = tmpDst[2].re + tmpDst[10].re;
    secondTmpDst[8].im = tmpDst[2].im + tmpDst[10].im;
    secondTmpDst[9].re = tmpDst[2].re - tmpDst[10].re;
    secondTmpDst[9].im = tmpDst[2].im - tmpDst[10].im;
    secondTmpDst[10].re = tmpDst[6].re + tmpDst[14].re;
    secondTmpDst[10].im = tmpDst[6].im + tmpDst[14].im;
    secondTmpDst[11].re = tmpDst[6].re - tmpDst[14].re;
    secondTmpDst[11].im = tmpDst[6].im - tmpDst[14].im;

    secondTmpDst[12].re = tmpDst[3].re + tmpDst[11].re;
    secondTmpDst[12].im = tmpDst[3].im + tmpDst[11].im;
    secondTmpDst[13].re = tmpDst[3].re - tmpDst[11].re;
    secondTmpDst[13].im = tmpDst[3].im - tmpDst[11].im;
    secondTmpDst[14].re = tmpDst[7].re + tmpDst[15].re;
    secondTmpDst[14].im = tmpDst[7].im + tmpDst[15].im;
    secondTmpDst[15].re = tmpDst[7].re - tmpDst[15].re;
    secondTmpDst[15].im = tmpDst[7].im - tmpDst[15].im;

    pDst[0].re = secondTmpDst[0].re + secondTmpDst[2].re;
    pDst[0].im = secondTmpDst[0].im + secondTmpDst[2].im;
    pDst[4].re = secondTmpDst[1].re + secondTmpDst[3].im;
    pDst[4].im = secondTmpDst[1].im - secondTmpDst[3].re;
    pDst[8].re = secondTmpDst[0].re - secondTmpDst[2].re;
    pDst[8].im = secondTmpDst[0].im - secondTmpDst[2].im;
    pDst[12].re = secondTmpDst[1].re - secondTmpDst[3].im;
    pDst[12].im = secondTmpDst[1].im + secondTmpDst[3].re;

    pDst[1].re = secondTmpDst[4].re + secondTmpDst[6].re;
    pDst[1].im = secondTmpDst[4].im + secondTmpDst[6].im;
    pDst[5].re = secondTmpDst[5].re + secondTmpDst[7].im;
    pDst[5].im = secondTmpDst[5].im - secondTmpDst[7].re;
    pDst[9].re = secondTmpDst[4].re - secondTmpDst[6].re;
    pDst[9].im = secondTmpDst[4].im - secondTmpDst[6].im;
    pDst[13].re = secondTmpDst[5].re - secondTmpDst[7].im;
    pDst[13].im = secondTmpDst[5].im + secondTmpDst[7].re;

    pDst[2].re = secondTmpDst[8].re + secondTmpDst[10].re;
    pDst[2].im = secondTmpDst[8].im + secondTmpDst[10].im;
    pDst[6].re = secondTmpDst[9].re + secondTmpDst[11].im;
    pDst[6].im = secondTmpDst[9].im - secondTmpDst[11].re;
    pDst[10].re = secondTmpDst[8].re - secondTmpDst[10].re;
    pDst[10].im = secondTmpDst[8].im - secondTmpDst[10].im;
    pDst[14].re = secondTmpDst[9].re - secondTmpDst[11].im;
    pDst[14].im = secondTmpDst[9].im + secondTmpDst[11].re;

    pDst[3].re = secondTmpDst[12].re + secondTmpDst[14].re;
    pDst[3].im = secondTmpDst[12].im + secondTmpDst[14].im;
    pDst[7].re = secondTmpDst[13].re + secondTmpDst[15].im;
    pDst[7].im = secondTmpDst[13].im - secondTmpDst[15].re;
    pDst[11].re = secondTmpDst[12].re - secondTmpDst[14].re;
    pDst[11].im = secondTmpDst[12].im - secondTmpDst[14].im;
    pDst[15].re = secondTmpDst[13].re - secondTmpDst[15].im;
    pDst[15].im = secondTmpDst[13].im + secondTmpDst[15].re;
    return;
}
