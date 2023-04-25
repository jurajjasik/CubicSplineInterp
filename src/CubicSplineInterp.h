#ifndef CubicSplineInterp_h
#define CubicSplineInterp_h

#include <Arduino.h>

#define CSI_MAX_TAB_POINTS 64

// Documentation: Numerical Recipes in C - The Art of Scientific Computing, 
// W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery, CAMBRIDGE UNIVERSITY PRESS, 1992,
// ISBN:978-0-521-43108-8

class CubicSplineInterp {
private:
	int _n; // number of tabulated points
	float _y2[CSI_MAX_TAB_POINTS]; // contains the second derivatives of the interpolating function at the tabulated points 
	const float* _x;
	const float* _y;
    
    int _klo;
    void _hunt(float x, int *jlo);
	
public:
	/**
	 *  \brief Brief
	 *  
	 *  \param [in] x[]
	 *  \param [in] y[]
	 *  \param [in] n Length of x[] and y[]
	 *  \return 
	 *  
	 *  \details Details
	 */
	void init(const float* x, const float* y, int n);
	
	float calc(float x);
    float calcHunt(float x);
};

#endif
