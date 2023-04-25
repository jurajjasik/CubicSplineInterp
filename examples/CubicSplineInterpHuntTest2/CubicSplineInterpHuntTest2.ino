#include <CubicSplineInterp.h>
#include <stdlib.h>

const float xt[] = {-1.0, 0.3, 2.4, 3.8};
const float yt[] = {1.0, 2.2, 0.2, 2.0};
size_t n = 4;

CubicSplineInterp spline = CubicSplineInterp();

void setup() {
	Serial.begin(115200);

	spline.init(xt, yt, n);
    
    for(int i=0; i<1000; ++i){
        float x = random(100000) / 100000.0 * 20.0 - 10.0; 
        Serial.print(x);
		Serial.print(" ");
		Serial.print( spline.calcHunt(x) );
		Serial.println();
    }
}

void loop() {

}