#include "lookup_table.h"
#include <math.h>


#define SIN W/4
#define INVPI 0.318309886183791
#define INV_TWOPI 0.159154943091896
//#define C3 0.0496389686056594 //1/((N+1)*BW*PI)
//#define N 512
#define C3 0.045391784126031
#define N 560

#define PI 3.14159265358979323846

extern float cosbuf[W];
extern float *sincbuf;

double cos_lut(double input){
	int index, trunc;
	double output;

	//Essentially, takes input MOD 2PI.
	input = input * INV_TWOPI;
	trunc = (int)input;
	input -= trunc;

	//Converts the input to the index of the corresponding value in cosbuf
	if(input >= 0)
		input *= W;
	else
		input = W - -W * input;

	/// May be better to round here.
	index = (int)input;


	if(index < 0)
		index += W;
	if(index >= W)
		index -= W;

	output = cosbuf[index];
	return output;
}



double sin_lut(double input){
	int index = 0, trunc = 0;
	double output = 0;

	input = input * INV_TWOPI; //W / (2 * PI)
	trunc = (int)input;
	input -= trunc;

	if(input >= 0)
		input *= W;
	else
		input = W - -W * input;

	//Was having problems adding in the sin shift before multiplying by W.
	index = (int)(input + 0.5) - SIN;

	if(index < 0)
		index += W;

	output = cosbuf[index];
	return output;
}



double sinc_lut(double input) {
	double output = 0;
	int index;

	//Input can be in between -N * BW and N * BW
	input = input * PI;

	input *= C3; // C3 is 1/((N)*BW*PI) //This results in a number in between -1 and 1
	input++; //Make the number between 0 and 2

	//A safety check, should never happen.
	if(input < 0)
		input = 0;

	input = 0.5 * W * input; //Make the number between 0 and W - 1

	index = (int)input;

	if(index >= W)
		index -= W;
	else if(index < 0)
		index += W;

	output = sincbuf[index];
	return output;
}


