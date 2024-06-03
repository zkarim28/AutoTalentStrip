//#ifndef MAYER_H
//#define MAYER_H
//
//#define REAL float
//
//void mayer_realfft(int n, REAL *real);
//void mayer_realifft(int n, REAL *real);
//
//#endif

#ifndef MAYER_H
#define MAYER_H

#define REAL float

#ifdef __cplusplus
extern "C" {
#endif

void mayer_realfft(int n, REAL *real);
void mayer_realifft(int n, REAL *real);

#ifdef __cplusplus
}
#endif

#endif // MAYER_H
