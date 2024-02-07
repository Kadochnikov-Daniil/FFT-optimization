#include "dft.h"

void dft32Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    //dft8 x dft4
    //DFT8 * 4
    cfloat32_t firstTmpDst[32];
    firstTmpDst[0].re = pSrc[0].re + pSrc[16].re;
    firstTmpDst[0].im = pSrc[0].im + pSrc[16].im;
    firstTmpDst[1].re = pSrc[0].re - pSrc[16].re;
    firstTmpDst[1].im = pSrc[0].im - pSrc[16].im;
    firstTmpDst[2].re = pSrc[4].re + pSrc[20].re;
    firstTmpDst[2].im = pSrc[4].im + pSrc[20].im;
    firstTmpDst[3].re = pSrc[4].re - pSrc[20].re;
    firstTmpDst[3].im = pSrc[4].im - pSrc[20].im;
    firstTmpDst[4].re = pSrc[8].re + pSrc[24].re;
    firstTmpDst[4].im = pSrc[8].im + pSrc[24].im;
    firstTmpDst[5].re = pSrc[8].im - pSrc[24].im;
    firstTmpDst[5].im = pSrc[24].re - pSrc[8].re;
    firstTmpDst[6].re = pSrc[12].re + pSrc[28].re;
    firstTmpDst[6].im = pSrc[12].im + pSrc[28].im;
    firstTmpDst[7].re = pSrc[28].re - pSrc[12].re;
    firstTmpDst[7].im = pSrc[28].im - pSrc[12].im;

    firstTmpDst[8].re = pSrc[1].re + pSrc[17].re;
    firstTmpDst[8].im = pSrc[1].im + pSrc[17].im;
    firstTmpDst[9].re = pSrc[1].re - pSrc[17].re;
    firstTmpDst[9].im = pSrc[1].im - pSrc[17].im;
    firstTmpDst[10].re = pSrc[5].re + pSrc[21].re;
    firstTmpDst[10].im = pSrc[5].im + pSrc[21].im;
    firstTmpDst[11].re = pSrc[5].re - pSrc[21].re;
    firstTmpDst[11].im = pSrc[5].im - pSrc[21].im;
    firstTmpDst[12].re = pSrc[9].re + pSrc[25].re;
    firstTmpDst[12].im = pSrc[9].im + pSrc[25].im;
    firstTmpDst[13].re = pSrc[9].im - pSrc[25].im;
    firstTmpDst[13].im = pSrc[25].re - pSrc[9].re;
    firstTmpDst[14].re = pSrc[13].re + pSrc[29].re;
    firstTmpDst[14].im = pSrc[13].im + pSrc[29].im;
    firstTmpDst[15].re = pSrc[29].re - pSrc[13].re;
    firstTmpDst[15].im = pSrc[29].im - pSrc[13].im;

    firstTmpDst[16].re = pSrc[2].re + pSrc[18].re;
    firstTmpDst[16].im = pSrc[2].im + pSrc[18].im;
    firstTmpDst[17].re = pSrc[2].re - pSrc[18].re;
    firstTmpDst[17].im = pSrc[2].im - pSrc[18].im;
    firstTmpDst[18].re = pSrc[6].re + pSrc[22].re;
    firstTmpDst[18].im = pSrc[6].im + pSrc[22].im;
    firstTmpDst[19].re = pSrc[6].re - pSrc[22].re;
    firstTmpDst[19].im = pSrc[6].im - pSrc[22].im;
    firstTmpDst[20].re = pSrc[10].re + pSrc[26].re;
    firstTmpDst[20].im = pSrc[10].im + pSrc[26].im;
    firstTmpDst[21].re = pSrc[10].im - pSrc[26].im;
    firstTmpDst[21].im = pSrc[26].re - pSrc[10].re;
    firstTmpDst[22].re = pSrc[14].re + pSrc[30].re;
    firstTmpDst[22].im = pSrc[14].im + pSrc[30].im;
    firstTmpDst[23].re = pSrc[30].re - pSrc[14].re;
    firstTmpDst[23].im = pSrc[30].im - pSrc[14].im;

    firstTmpDst[24].re = pSrc[3].re + pSrc[19].re;
    firstTmpDst[24].im = pSrc[3].im + pSrc[19].im;
    firstTmpDst[25].re = pSrc[3].re - pSrc[19].re;
    firstTmpDst[25].im = pSrc[3].im - pSrc[19].im;
    firstTmpDst[26].re = pSrc[7].re + pSrc[23].re;
    firstTmpDst[26].im = pSrc[7].im + pSrc[23].im;
    firstTmpDst[27].re = pSrc[7].re - pSrc[23].re;
    firstTmpDst[27].im = pSrc[7].im - pSrc[23].im;
    firstTmpDst[28].re = pSrc[11].re + pSrc[27].re;
    firstTmpDst[28].im = pSrc[11].im + pSrc[27].im;
    firstTmpDst[29].re = pSrc[11].im - pSrc[27].im;
    firstTmpDst[29].im = pSrc[27].re - pSrc[11].re;
    firstTmpDst[30].re = pSrc[15].re + pSrc[31].re;
    firstTmpDst[30].im = pSrc[15].im + pSrc[31].im;
    firstTmpDst[31].re = pSrc[31].re - pSrc[15].re;
    firstTmpDst[31].im = pSrc[31].im - pSrc[15].im;

    //coefficient multiplication for DFT8
    float pi_4 = 0.707106781186548;
    float im = (firstTmpDst[3].im - firstTmpDst[3].re) * pi_4;
    firstTmpDst[3].re = (firstTmpDst[3].re + firstTmpDst[3].im) * pi_4;
    firstTmpDst[3].im = im;
    im = (firstTmpDst[11].im - firstTmpDst[11].re) * pi_4;
    firstTmpDst[11].re = (firstTmpDst[11].re + firstTmpDst[11].im) * pi_4;
    firstTmpDst[11].im = im;
    im = (firstTmpDst[19].im - firstTmpDst[19].re) * pi_4;
    firstTmpDst[19].re = (firstTmpDst[19].re + firstTmpDst[19].im) * pi_4;
    firstTmpDst[19].im = im;
    im = (firstTmpDst[27].im - firstTmpDst[27].re) * pi_4;
    firstTmpDst[27].re = (firstTmpDst[27].re + firstTmpDst[27].im) * pi_4;
    firstTmpDst[27].im = im;

    im = (firstTmpDst[7].im + firstTmpDst[7].re) * pi_4;
    firstTmpDst[7].re = (firstTmpDst[7].re - firstTmpDst[7].im) * pi_4;
    firstTmpDst[7].im = im;
    im = (firstTmpDst[15].im + firstTmpDst[15].re) * pi_4;
    firstTmpDst[15].re = (firstTmpDst[15].re - firstTmpDst[15].im) * pi_4;
    firstTmpDst[15].im = im;
    im = (firstTmpDst[23].im + firstTmpDst[23].re) * pi_4;
    firstTmpDst[23].re = (firstTmpDst[23].re - firstTmpDst[23].im) * pi_4;
    firstTmpDst[23].im = im;
    im = (firstTmpDst[31].im + firstTmpDst[31].re) * pi_4;
    firstTmpDst[31].re = (firstTmpDst[31].re - firstTmpDst[31].im) * pi_4;
    firstTmpDst[31].im = im;

    //DFT4 for DFT8
    cfloat32_t secondTmpDst[32];
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

    secondTmpDst[8].re = firstTmpDst[8].re + firstTmpDst[12].re;
    secondTmpDst[8].im = firstTmpDst[8].im + firstTmpDst[12].im;
    secondTmpDst[9].re = firstTmpDst[8].re - firstTmpDst[12].re;
    secondTmpDst[9].im = firstTmpDst[8].im - firstTmpDst[12].im;
    secondTmpDst[10].re = firstTmpDst[10].re + firstTmpDst[14].re;
    secondTmpDst[10].im = firstTmpDst[10].im + firstTmpDst[14].im;
    secondTmpDst[11].re = firstTmpDst[10].re - firstTmpDst[14].re;
    secondTmpDst[11].im = firstTmpDst[10].im - firstTmpDst[14].im;
    secondTmpDst[12].re = firstTmpDst[9].re + firstTmpDst[13].re;
    secondTmpDst[12].im = firstTmpDst[9].im + firstTmpDst[13].im;
    secondTmpDst[13].re = firstTmpDst[9].re - firstTmpDst[13].re;
    secondTmpDst[13].im = firstTmpDst[9].im - firstTmpDst[13].im;
    secondTmpDst[14].re = firstTmpDst[11].re + firstTmpDst[15].re;
    secondTmpDst[14].im = firstTmpDst[11].im + firstTmpDst[15].im;
    secondTmpDst[15].re = firstTmpDst[11].re - firstTmpDst[15].re;
    secondTmpDst[15].im = firstTmpDst[11].im - firstTmpDst[15].im;

    secondTmpDst[16].re = firstTmpDst[16].re + firstTmpDst[20].re;
    secondTmpDst[16].im = firstTmpDst[16].im + firstTmpDst[20].im;
    secondTmpDst[17].re = firstTmpDst[16].re - firstTmpDst[20].re;
    secondTmpDst[17].im = firstTmpDst[16].im - firstTmpDst[20].im;
    secondTmpDst[18].re = firstTmpDst[18].re + firstTmpDst[22].re;
    secondTmpDst[18].im = firstTmpDst[18].im + firstTmpDst[22].im;
    secondTmpDst[19].re = firstTmpDst[18].re - firstTmpDst[22].re;
    secondTmpDst[19].im = firstTmpDst[18].im - firstTmpDst[22].im;
    secondTmpDst[20].re = firstTmpDst[17].re + firstTmpDst[21].re;
    secondTmpDst[20].im = firstTmpDst[17].im + firstTmpDst[21].im;
    secondTmpDst[21].re = firstTmpDst[17].re - firstTmpDst[21].re;
    secondTmpDst[21].im = firstTmpDst[17].im - firstTmpDst[21].im;
    secondTmpDst[22].re = firstTmpDst[19].re + firstTmpDst[23].re;
    secondTmpDst[22].im = firstTmpDst[19].im + firstTmpDst[23].im;
    secondTmpDst[23].re = firstTmpDst[19].re - firstTmpDst[23].re;
    secondTmpDst[23].im = firstTmpDst[19].im - firstTmpDst[23].im;

    secondTmpDst[24].re = firstTmpDst[24].re + firstTmpDst[28].re;
    secondTmpDst[24].im = firstTmpDst[24].im + firstTmpDst[28].im;
    secondTmpDst[25].re = firstTmpDst[24].re - firstTmpDst[28].re;
    secondTmpDst[25].im = firstTmpDst[24].im - firstTmpDst[28].im;
    secondTmpDst[26].re = firstTmpDst[26].re + firstTmpDst[30].re;
    secondTmpDst[26].im = firstTmpDst[26].im + firstTmpDst[30].im;
    secondTmpDst[27].re = firstTmpDst[30].re - firstTmpDst[26].re;
    secondTmpDst[27].im = firstTmpDst[26].im - firstTmpDst[30].im;
    secondTmpDst[28].re = firstTmpDst[25].re + firstTmpDst[29].re;
    secondTmpDst[28].im = firstTmpDst[25].im + firstTmpDst[29].im;
    secondTmpDst[29].re = firstTmpDst[25].re - firstTmpDst[29].re;
    secondTmpDst[29].im = firstTmpDst[25].im - firstTmpDst[29].im;
    secondTmpDst[30].re = firstTmpDst[27].re + firstTmpDst[31].re;
    secondTmpDst[30].im = firstTmpDst[27].im + firstTmpDst[31].im;
    secondTmpDst[31].re = firstTmpDst[31].re - firstTmpDst[27].re;
    secondTmpDst[31].im = firstTmpDst[31].im - firstTmpDst[27].im;

    cfloat32_t tmpDst[32];
    tmpDst[0].re = secondTmpDst[0].re + secondTmpDst[2].re;
    tmpDst[0].im = secondTmpDst[0].im + secondTmpDst[2].im;
    tmpDst[2].re = secondTmpDst[1].re + secondTmpDst[3].im;
    tmpDst[2].im = secondTmpDst[1].im - secondTmpDst[3].re;
    tmpDst[4].re = secondTmpDst[0].re - secondTmpDst[2].re;
    tmpDst[4].im = secondTmpDst[0].im - secondTmpDst[2].im;
    tmpDst[6].re = secondTmpDst[1].re - secondTmpDst[3].im;
    tmpDst[6].im = secondTmpDst[1].im + secondTmpDst[3].re;
    tmpDst[1].re = secondTmpDst[4].re + secondTmpDst[6].re;
    tmpDst[1].im = secondTmpDst[4].im + secondTmpDst[6].im;
    tmpDst[3].re = secondTmpDst[5].re + secondTmpDst[7].im;
    tmpDst[3].im = secondTmpDst[5].im - secondTmpDst[7].re;
    tmpDst[5].re = secondTmpDst[4].re - secondTmpDst[6].re;
    tmpDst[5].im = secondTmpDst[4].im - secondTmpDst[6].im;
    tmpDst[7].re = secondTmpDst[5].re - secondTmpDst[7].im;
    tmpDst[7].im = secondTmpDst[5].im + secondTmpDst[7].re;

    tmpDst[8].re = secondTmpDst[8].re + secondTmpDst[10].re;
    tmpDst[8].im = secondTmpDst[8].im + secondTmpDst[10].im;
    tmpDst[10].re = secondTmpDst[9].re + secondTmpDst[11].im;
    tmpDst[10].im = secondTmpDst[9].im - secondTmpDst[11].re;
    tmpDst[12].re = secondTmpDst[8].re - secondTmpDst[10].re;
    tmpDst[12].im = secondTmpDst[8].im - secondTmpDst[10].im;
    tmpDst[14].re = secondTmpDst[9].re - secondTmpDst[11].im;
    tmpDst[14].im = secondTmpDst[9].im + secondTmpDst[11].re;
    tmpDst[9].re = secondTmpDst[12].re + secondTmpDst[14].re;
    tmpDst[9].im = secondTmpDst[12].im + secondTmpDst[14].im;
    tmpDst[11].re = secondTmpDst[13].re + secondTmpDst[15].im;
    tmpDst[11].im = secondTmpDst[13].im - secondTmpDst[15].re;
    tmpDst[13].re = secondTmpDst[12].re - secondTmpDst[14].re;
    tmpDst[13].im = secondTmpDst[12].im - secondTmpDst[14].im;
    tmpDst[15].re = secondTmpDst[13].re - secondTmpDst[15].im;
    tmpDst[15].im = secondTmpDst[13].im + secondTmpDst[15].re;

    tmpDst[16].re = secondTmpDst[16].re + secondTmpDst[18].re;
    tmpDst[16].im = secondTmpDst[16].im + secondTmpDst[18].im;
    tmpDst[18].re = secondTmpDst[17].re + secondTmpDst[19].im;
    tmpDst[18].im = secondTmpDst[17].im - secondTmpDst[19].re;
    tmpDst[20].re = secondTmpDst[16].im - secondTmpDst[18].im;
    tmpDst[20].im = secondTmpDst[18].re - secondTmpDst[16].re;
    tmpDst[22].re = secondTmpDst[19].im - secondTmpDst[17].re;
    tmpDst[22].im = secondTmpDst[17].im + secondTmpDst[19].re;
    tmpDst[17].re = secondTmpDst[20].re + secondTmpDst[22].re;
    tmpDst[17].im = secondTmpDst[20].im + secondTmpDst[22].im;
    tmpDst[19].re = secondTmpDst[21].re + secondTmpDst[23].im;
    tmpDst[19].im = secondTmpDst[21].im - secondTmpDst[23].re;
    tmpDst[21].re = secondTmpDst[22].re - secondTmpDst[20].re;
    tmpDst[21].im = secondTmpDst[22].im - secondTmpDst[20].im;
    tmpDst[23].re = secondTmpDst[23].im - secondTmpDst[21].re;
    tmpDst[23].im = secondTmpDst[21].im + secondTmpDst[23].re;

    tmpDst[24].re = secondTmpDst[24].re + secondTmpDst[26].re;
    tmpDst[24].im = secondTmpDst[24].im + secondTmpDst[26].im;
    tmpDst[26].re = secondTmpDst[25].re + secondTmpDst[27].im;
    tmpDst[26].im = secondTmpDst[25].im + secondTmpDst[27].re;
    tmpDst[28].re = secondTmpDst[26].re - secondTmpDst[24].re;
    tmpDst[28].im = secondTmpDst[24].im - secondTmpDst[26].im;
    tmpDst[30].re = secondTmpDst[25].re - secondTmpDst[27].im;
    tmpDst[30].im = secondTmpDst[27].re - secondTmpDst[25].im;
    tmpDst[25].re = secondTmpDst[28].re + secondTmpDst[30].re;
    tmpDst[25].im = secondTmpDst[28].im + secondTmpDst[30].im;
    tmpDst[27].re = secondTmpDst[31].im - secondTmpDst[29].re;
    tmpDst[27].im = secondTmpDst[29].im + secondTmpDst[31].re;
    tmpDst[29].re = secondTmpDst[30].re - secondTmpDst[28].re;
    tmpDst[29].im = secondTmpDst[28].im - secondTmpDst[30].im;
    tmpDst[31].re = secondTmpDst[29].re + secondTmpDst[31].im;
    tmpDst[31].im = secondTmpDst[31].re - secondTmpDst[29].im;
    
    //coefficient multiplication for dft32
    float pi_16 = 0.980785250663757;
    float pi_8 = 0.923879504203796;
    float pi3_16 = 0.831469595432281;
    //float pi_4 = 0.707106769084930; //defined in upper DFT8
    float pi5_16 = 0.555570244789124;
    float pi3_8 = 0.382683426141739;
    float pi7_16 = 0.195090323686600;

    im = tmpDst[9].im * pi_16 - tmpDst[9].re * pi7_16;
    tmpDst[9].re = tmpDst[9].re * pi_16 + tmpDst[9].im * pi7_16;
    tmpDst[9].im = im;

    im = tmpDst[10].im * pi_8 - tmpDst[10].re * pi3_8;
    tmpDst[10].re = tmpDst[10].re * pi_8 + tmpDst[10].im * pi3_8;
    tmpDst[10].im = im;

    im = tmpDst[11].im * pi3_16 - tmpDst[11].re * pi5_16;
    tmpDst[11].re = tmpDst[11].re * pi3_16 + tmpDst[11].im * pi5_16;
    tmpDst[11].im = im;

    im = (tmpDst[12].im - tmpDst[12].re) * pi_4;
    tmpDst[12].re = (tmpDst[12].re + tmpDst[12].im) * pi_4;
    tmpDst[12].im = im;

    im = tmpDst[13].im * pi5_16 - tmpDst[13].re * pi3_16;
    tmpDst[13].re = tmpDst[13].re * pi5_16 + tmpDst[13].im * pi3_16;
    tmpDst[13].im = im;

    im = tmpDst[14].im * pi3_8 - tmpDst[14].re * pi_8;
    tmpDst[14].re = tmpDst[14].re * pi3_8 + tmpDst[14].im * pi_8;
    tmpDst[14].im = im;

    im = tmpDst[15].im * pi7_16 - tmpDst[15].re * pi_16;
    tmpDst[15].re = tmpDst[15].re * pi7_16 + tmpDst[15].im * pi_16;
    tmpDst[15].im = im;

    im = tmpDst[17].im * pi_8 - tmpDst[17].re * pi3_8;
    tmpDst[17].re = tmpDst[17].re * pi_8 + tmpDst[17].im * pi3_8;
    tmpDst[17].im = im;

    im = (tmpDst[18].im - tmpDst[18].re) * pi_4;
    tmpDst[18].re = (tmpDst[18].re + tmpDst[18].im) * pi_4;
    tmpDst[18].im = im;

    im = tmpDst[19].im * pi3_8 - tmpDst[19].re * pi_8;
    tmpDst[19].re = tmpDst[19].re * pi3_8 + tmpDst[19].im * pi_8;
    tmpDst[19].im = im;

    im = tmpDst[21].im * pi3_8 + tmpDst[21].re * pi_8;
    tmpDst[21].re = tmpDst[21].re * pi3_8 - tmpDst[21].im * pi_8;
    tmpDst[21].im = im;

    im = (tmpDst[22].re - tmpDst[22].im) * pi_4;
    tmpDst[22].re = (tmpDst[22].im + tmpDst[22].re) * pi_4;
    tmpDst[22].im = im;

    im = tmpDst[23].re * pi3_8 - tmpDst[23].im * pi_8;
    tmpDst[23].re = tmpDst[23].im * pi3_8 + tmpDst[23].re * pi_8;
    tmpDst[23].im = im;

    im = tmpDst[25].im * pi3_16 - tmpDst[25].re * pi5_16;
    tmpDst[25].re = tmpDst[25].re * pi3_16 + tmpDst[25].im * pi5_16;
    tmpDst[25].im = im;

    im = tmpDst[26].im * pi3_8 - tmpDst[26].re * pi_8;
    tmpDst[26].re = tmpDst[26].re * pi3_8 + tmpDst[26].im * pi_8;
    tmpDst[26].im = im;

    im = tmpDst[27].re * pi_16 - tmpDst[27].im * pi7_16;
    tmpDst[27].re = tmpDst[27].im * pi_16 + tmpDst[27].re * pi7_16;
    tmpDst[27].im = im;

    im = (tmpDst[28].re - tmpDst[28].im) * pi_4;
    tmpDst[28].re = (tmpDst[28].re + tmpDst[28].im) * pi_4;
    tmpDst[28].im = im;

    im = tmpDst[29].re * pi7_16 - tmpDst[29].im * pi_16;
    tmpDst[29].re = tmpDst[29].re * pi_16 + tmpDst[29].im * pi7_16;
    tmpDst[29].im = im;

    im = tmpDst[30].im * pi_8 + tmpDst[30].re * pi3_8;
    tmpDst[30].re = tmpDst[30].im * pi3_8 - tmpDst[30].re * pi_8;
    tmpDst[30].im = im;

    im = tmpDst[31].im * pi5_16 + tmpDst[31].re * pi3_16;
    tmpDst[31].re = tmpDst[31].im * pi3_16 - tmpDst[31].re * pi5_16;
    tmpDst[31].im = im;

    //DFT4 * 8
    secondTmpDst[0].re = tmpDst[0].re + tmpDst[16].re;
    secondTmpDst[0].im = tmpDst[0].im + tmpDst[16].im;
    secondTmpDst[1].re = tmpDst[0].re - tmpDst[16].re;
    secondTmpDst[1].im = tmpDst[0].im - tmpDst[16].im;
    secondTmpDst[2].re = tmpDst[8].re + tmpDst[24].re;
    secondTmpDst[2].im = tmpDst[8].im + tmpDst[24].im;
    secondTmpDst[3].re = tmpDst[8].re - tmpDst[24].re;
    secondTmpDst[3].im = tmpDst[8].im - tmpDst[24].im;

    secondTmpDst[4].re = tmpDst[1].re + tmpDst[17].re;
    secondTmpDst[4].im = tmpDst[1].im + tmpDst[17].im;
    secondTmpDst[5].re = tmpDst[1].re - tmpDst[17].re;
    secondTmpDst[5].im = tmpDst[1].im - tmpDst[17].im;
    secondTmpDst[6].re = tmpDst[9].re + tmpDst[25].re;
    secondTmpDst[6].im = tmpDst[9].im + tmpDst[25].im;
    secondTmpDst[7].re = tmpDst[9].re - tmpDst[25].re;
    secondTmpDst[7].im = tmpDst[9].im - tmpDst[25].im;

    secondTmpDst[8].re = tmpDst[2].re + tmpDst[18].re;
    secondTmpDst[8].im = tmpDst[2].im + tmpDst[18].im;
    secondTmpDst[9].re = tmpDst[2].re - tmpDst[18].re;
    secondTmpDst[9].im = tmpDst[2].im - tmpDst[18].im;
    secondTmpDst[10].re = tmpDst[10].re + tmpDst[26].re;
    secondTmpDst[10].im = tmpDst[10].im + tmpDst[26].im;
    secondTmpDst[11].re = tmpDst[10].re - tmpDst[26].re;
    secondTmpDst[11].im = tmpDst[10].im - tmpDst[26].im;

    secondTmpDst[12].re = tmpDst[3].re + tmpDst[19].re;
    secondTmpDst[12].im = tmpDst[3].im + tmpDst[19].im;
    secondTmpDst[13].re = tmpDst[3].re - tmpDst[19].re;
    secondTmpDst[13].im = tmpDst[3].im - tmpDst[19].im;
    secondTmpDst[14].re = tmpDst[11].re + tmpDst[27].re;
    secondTmpDst[14].im = tmpDst[11].im + tmpDst[27].im;
    secondTmpDst[15].re = tmpDst[11].re - tmpDst[27].re;
    secondTmpDst[15].im = tmpDst[11].im - tmpDst[27].im;

    secondTmpDst[16].re = tmpDst[4].re + tmpDst[20].re;
    secondTmpDst[16].im = tmpDst[4].im + tmpDst[20].im;
    secondTmpDst[17].re = tmpDst[4].re - tmpDst[20].re;
    secondTmpDst[17].im = tmpDst[4].im - tmpDst[20].im;
    secondTmpDst[18].re = tmpDst[12].re + tmpDst[28].re;
    secondTmpDst[18].im = tmpDst[12].im + tmpDst[28].im;
    secondTmpDst[19].re = tmpDst[12].re - tmpDst[28].re;
    secondTmpDst[19].im = tmpDst[12].im - tmpDst[28].im;

    secondTmpDst[20].re = tmpDst[5].re + tmpDst[21].re;
    secondTmpDst[20].im = tmpDst[5].im + tmpDst[21].im;
    secondTmpDst[21].re = tmpDst[5].re - tmpDst[21].re;
    secondTmpDst[21].im = tmpDst[5].im - tmpDst[21].im;
    secondTmpDst[22].re = tmpDst[13].re + tmpDst[29].re;
    secondTmpDst[22].im = tmpDst[13].im + tmpDst[29].im;
    secondTmpDst[23].re = tmpDst[13].re - tmpDst[29].re;
    secondTmpDst[23].im = tmpDst[13].im - tmpDst[29].im;

    secondTmpDst[24].re = tmpDst[6].re + tmpDst[22].re;
    secondTmpDst[24].im = tmpDst[6].im + tmpDst[22].im;
    secondTmpDst[25].re = tmpDst[6].re - tmpDst[22].re;
    secondTmpDst[25].im = tmpDst[6].im - tmpDst[22].im;
    secondTmpDst[26].re = tmpDst[14].re + tmpDst[30].re;
    secondTmpDst[26].im = tmpDst[14].im + tmpDst[30].im;
    secondTmpDst[27].re = tmpDst[14].re - tmpDst[30].re;
    secondTmpDst[27].im = tmpDst[14].im - tmpDst[30].im;

    secondTmpDst[28].re = tmpDst[7].re + tmpDst[23].re;
    secondTmpDst[28].im = tmpDst[7].im + tmpDst[23].im;
    secondTmpDst[29].re = tmpDst[7].re - tmpDst[23].re;
    secondTmpDst[29].im = tmpDst[7].im - tmpDst[23].im;
    secondTmpDst[30].re = tmpDst[15].re + tmpDst[31].re;
    secondTmpDst[30].im = tmpDst[15].im + tmpDst[31].im;
    secondTmpDst[31].re = tmpDst[15].re - tmpDst[31].re;
    secondTmpDst[31].im = tmpDst[15].im - tmpDst[31].im;

    pDst[0].re = secondTmpDst[0].re + secondTmpDst[2].re;
    pDst[0].im = secondTmpDst[0].im + secondTmpDst[2].im;
    pDst[8].re = secondTmpDst[1].re + secondTmpDst[3].im;
    pDst[8].im = secondTmpDst[1].im - secondTmpDst[3].re;
    pDst[16].re = secondTmpDst[0].re - secondTmpDst[2].re;
    pDst[16].im = secondTmpDst[0].im - secondTmpDst[2].im;
    pDst[24].re = secondTmpDst[1].re - secondTmpDst[3].im;
    pDst[24].im = secondTmpDst[1].im + secondTmpDst[3].re;

    pDst[1].re = secondTmpDst[4].re + secondTmpDst[6].re;
    pDst[1].im = secondTmpDst[4].im + secondTmpDst[6].im;
    pDst[9].re = secondTmpDst[5].re + secondTmpDst[7].im;
    pDst[9].im = secondTmpDst[5].im - secondTmpDst[7].re;
    pDst[17].re = secondTmpDst[4].re - secondTmpDst[6].re;
    pDst[17].im = secondTmpDst[4].im - secondTmpDst[6].im;
    pDst[25].re = secondTmpDst[5].re - secondTmpDst[7].im;
    pDst[25].im = secondTmpDst[5].im + secondTmpDst[7].re;
    
    pDst[2].re = secondTmpDst[8].re + secondTmpDst[10].re;
    pDst[2].im = secondTmpDst[8].im + secondTmpDst[10].im;
    pDst[10].re = secondTmpDst[9].re + secondTmpDst[11].im;
    pDst[10].im = secondTmpDst[9].im - secondTmpDst[11].re;
    pDst[18].re = secondTmpDst[8].re - secondTmpDst[10].re;
    pDst[18].im = secondTmpDst[8].im - secondTmpDst[10].im;
    pDst[26].re = secondTmpDst[9].re - secondTmpDst[11].im;
    pDst[26].im = secondTmpDst[9].im + secondTmpDst[11].re;

    pDst[3].re = secondTmpDst[12].re + secondTmpDst[14].re;
    pDst[3].im = secondTmpDst[12].im + secondTmpDst[14].im;
    pDst[11].re = secondTmpDst[13].re + secondTmpDst[15].im;
    pDst[11].im = secondTmpDst[13].im - secondTmpDst[15].re;
    pDst[19].re = secondTmpDst[12].re - secondTmpDst[14].re;
    pDst[19].im = secondTmpDst[12].im - secondTmpDst[14].im;
    pDst[27].re = secondTmpDst[13].re - secondTmpDst[15].im;
    pDst[27].im = secondTmpDst[13].im + secondTmpDst[15].re;

    pDst[4].re = secondTmpDst[16].re + secondTmpDst[18].re;
    pDst[4].im = secondTmpDst[16].im + secondTmpDst[18].im;
    pDst[12].re = secondTmpDst[17].re + secondTmpDst[19].im;
    pDst[12].im = secondTmpDst[17].im - secondTmpDst[19].re;
    pDst[20].re = secondTmpDst[16].re - secondTmpDst[18].re;
    pDst[20].im = secondTmpDst[16].im - secondTmpDst[18].im;
    pDst[28].re = secondTmpDst[17].re - secondTmpDst[19].im;
    pDst[28].im = secondTmpDst[17].im + secondTmpDst[19].re;

    pDst[5].re = secondTmpDst[20].re + secondTmpDst[22].re;
    pDst[5].im = secondTmpDst[20].im + secondTmpDst[22].im;
    pDst[13].re = secondTmpDst[21].re + secondTmpDst[23].im;
    pDst[13].im = secondTmpDst[21].im - secondTmpDst[23].re;
    pDst[21].re = secondTmpDst[20].re - secondTmpDst[22].re;
    pDst[21].im = secondTmpDst[20].im - secondTmpDst[22].im;
    pDst[29].re = secondTmpDst[21].re - secondTmpDst[23].im;
    pDst[29].im = secondTmpDst[21].im + secondTmpDst[23].re;

    pDst[6].re = secondTmpDst[24].re + secondTmpDst[26].re;
    pDst[6].im = secondTmpDst[24].im + secondTmpDst[26].im;
    pDst[14].re = secondTmpDst[25].re + secondTmpDst[27].im;
    pDst[14].im = secondTmpDst[25].im - secondTmpDst[27].re;
    pDst[22].re = secondTmpDst[24].re - secondTmpDst[26].re;
    pDst[22].im = secondTmpDst[24].im - secondTmpDst[26].im;
    pDst[30].re = secondTmpDst[25].re - secondTmpDst[27].im;
    pDst[30].im = secondTmpDst[25].im + secondTmpDst[27].re;

    pDst[7].re = secondTmpDst[28].re + secondTmpDst[30].re;
    pDst[7].im = secondTmpDst[28].im + secondTmpDst[30].im;
    pDst[15].re = secondTmpDst[29].re + secondTmpDst[31].im;
    pDst[15].im = secondTmpDst[29].im - secondTmpDst[31].re;
    pDst[23].re = secondTmpDst[28].re - secondTmpDst[30].re;
    pDst[23].im = secondTmpDst[28].im - secondTmpDst[30].im;
    pDst[31].re = secondTmpDst[29].re - secondTmpDst[31].im;
    pDst[31].im = secondTmpDst[29].im + secondTmpDst[31].re;
    return;
}
