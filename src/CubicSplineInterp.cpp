#include "CubicSplineInterp.h"

#define CSI_TRACE(x)
//#define CSI_TRACE(x) x

void CubicSplineInterp::init(const float* x, const float* y, int n){
    while (_inCalc) {}

    _notValid = true;

    _x = x;
    _y = y;
    _n = n > CSI_MAX_TAB_POINTS ? CSI_MAX_TAB_POINTS : n;

    float u[_n];

    _y2[0] = u[0] = 0.0; // natural spline

    for (int i = 1; i < _n-1; ++i) {
        float sig = ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] );
        float p = sig * _y2[i-1] + 2.0;
        _y2[i] = ( sig - 1.0 ) / p;
        u[i] = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] )
            - ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
        u[i] = ( 6.0 * u[i] / ( x[i+1] - x[i-1] ) - sig * u[i-1] ) / p;
    }

    _y2[_n-1] = 0.0;

    for ( int i = _n-2; i > 0; --i) {
        _y2[i] = _y2[i] * _y2[i+1] + u[i];
    }

    _klo = 0;

    _notValid = false;
}

float CubicSplineInterp::calc(float x) const
{
    if (_notValid) return 0.0;

    int klo = 0;
    int khi = _n - 1;
    while ( khi - klo > 1 ) {
        int k = ( khi + klo ) >> 1;
        if ( _x[k] > x ) khi = k;
        else klo = k;
    }

    float h = _x[khi] - _x[klo];
    if ( !( h > 0.0 || h < 0.0 ) )
        return 0.0;
    float a = ( _x[khi] - x ) / h;
    float b = ( x - _x[klo] ) / h;
    return a * _y[klo] + b * _y[khi] +
        ( ( a * a * a - a ) * _y2[klo] + ( b * b * b - b ) * _y2[khi] ) * ( h * h ) / 6.0;
}

float CubicSplineInterp::calcHunt(float x)
{
    if (_notValid) return 0.0;

    _hunt(x, &_klo);

    int klo = _klo;

    if(klo >= _n - 1) klo = _n - 2;
    if(klo < 0 ) klo = 0;

    CSI_TRACE( printf("calcHunt => klo: %d\r\n", klo ); )

    float h = _x[klo + 1] - _x[klo];

    if ( !( h > 0.0 || h < 0.0 ) )
        return 0.0;
    float a = ( _x[klo + 1] - x ) / h;
    float b = ( x - _x[klo] ) / h;
    return a * _y[klo] + b * _y[klo + 1] +
        ( ( a * a * a - a ) * _y2[klo] + ( b * b * b - b ) * _y2[klo + 1] ) * ( h * h ) / 6.0;
}

void CubicSplineInterp::_hunt(float x, int *jlo){
    _inCalc = true;

    int jhi;
    int ascnd = ( _x[_n - 1] >= _x[0] ); //True if ascending order of table, false otherwise.

    if(ascnd){
        if (x > _x[_n - 1]) {
            *jlo = _n;
            return;
        }
        if (x < _x[0]) {
            *jlo = -1;
            return;
        }
    } else {
        if (x < _x[_n - 1]) {
            *jlo = _n;
            return;
        }
        if (x > _x[0]) {
            *jlo = -1;
            return;
        }
    }

    if (*jlo < 0 || *jlo >= _n) { //Input guess not useful. Go immediately to bisection
        *jlo = -1;
        jhi = _n - 1;
    } else {
        unsigned int inc = 1; // Set the hunting increment.
        if ( x >= _x[*jlo] == ascnd ) { // Hunt up:
            CSI_TRACE( printf("_hunt => hunt up\r\n"); )
            if (*jlo == _n - 1)
                return;
            jhi = (*jlo) + 1;
            while ( x >= _x[jhi] == ascnd ) { // Not done hunting,
                *jlo = jhi;
                inc <<= 1; // so double the increment
                jhi = (*jlo) + inc;
                if (jhi > _n - 1) { // Done hunting, since off end of table.
                    jhi = _n;
                    break;
                }
            } // Done hunting, value bracketed.
        } else { // Hunt down:
            CSI_TRACE( printf("_hunt => hunt down\r\n"); )
            if (*jlo == 0) {
                return;
            }
            jhi = (*jlo)--;
            while (x < _x[*jlo] == ascnd ) { // Not done hunting,
                jhi = (*jlo);
                inc <<= 1; // so double the increment
                if (inc >= jhi) { // Done hunting, since off end of table.
                    *jlo = 0;
                    break;
                }
                else *jlo = jhi - inc;
            } //and try again.
        } //Done hunting, value bracketed.
    } // Hunt is done, so begin the final bisection phase:
    CSI_TRACE( printf("_hunt => done hunting up, jlo:%d, jhi: %d\r\n", (*jlo), jhi); )
    while ( jhi - (*jlo) != 1) {
        int jm = ( jhi + (*jlo) ) >> 1;
        if (x >= _x[jm] == ascnd)
            *jlo = jm;
        else
            jhi = jm;
    }

    _inCalc = false;
}


