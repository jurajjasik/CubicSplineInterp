#include <CubicSplineInterp.h>

const float xt[] = {-1.0, 0.3, 2.4, 3.8};
const float yt[] = {1.0, 2.2, 0.2, 2.0};
size_t n = 4;

CubicSplineInterp spline = CubicSplineInterp();

void setup() {
	Serial.begin(115200);

	spline.init(xt, yt, n);
	
	for(float x = -10.0; x < 10.0; x += .1 ){
		Serial.print(x);
		Serial.print(" ");
		Serial.print( spline.calc(x) );
		Serial.println();
	}
}

void loop() {

}