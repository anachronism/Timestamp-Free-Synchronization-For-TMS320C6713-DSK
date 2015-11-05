/*************************************************************************
 *  Timestamp free synchronization
 *  Master node code
 *  DRB Apr 10, 2015
 *  Uses time reversal ideas developed on Apr 9, 2015
 *  States:
 *  0 : searching
 *  1 : recording
 *  2 : waiting for clock tick
 *  3 : waiting to play back response
 *  4 : playing back response
 *
 *  Updated October 25, 2015 to add a high pass filter on the input
 *  to attenuate signals at 1kHz carrier
 *
 *************************************************************************/

#define CHIP_6713 1

// clock pulse carrier kHz
// until October 2015, this was always 2
#define CARRIERKHZ 2

// which channel to output clock pulses
// until October 2015, this was always the left channel (0=right, 1=left)
#define CLOCKOUTPUTCHANNEL 1

// clock and sync pulse amplitude (should be no larger than 16383 if pulses are superimposed)
#define CLOCKAMPLITUDE 32000

// high pass filter *negative* group delay (HPF has group delay of 10)  // -14
#define HPFGROUPDELAY 0

// length of searching window in samples
#define M 40

// threshold value for searching window
#define T1 100000

// sinc pulse normalized bandwidth
#define BW 0.0125

// 2*N+1 is the number of samples in the sinc function
#define N 560 // was 512

// virtual clock counter
#define L 4096

// set this to one to pass through the left channel
// (useful for debugging master node functionality)
// or zero to not pass through (should normally be zero)
#define PASSTHROUGH 0

// define PI and INVPI
#define PI 3.14159265358979323846
#define INVPI 0.318309886183791

#include <stdio.h>
#include <c6x.h>
#include <csl.h>
#include <csl_mcbsp.h>
#include <csl_irq.h>
#include <math.h>

#include "dsk6713.h"
#include "dsk6713_aic23.h"
#include "dsk6713_led.h"

// ------------------------------------------
// start of variables
// ------------------------------------------
float buf[M];       // search buffer
float fbuf[M];      // filtered search buffer
float mfc[M];		// in-phase correlation buffer
float mfs[M];       // quadrature correlation buffer
float corr_max, corr_max_s, corr_max_c; // correlation variables
float corr_c[2*M];
float corr_s[2*M];
float s[2*M];
short corr_max_lag;
short bufindex = 0;
short i,j;
float zc, zs, z;
double t,x,y;
short wait_count = 0;
short clockbuf[L];  // clock buffer (right channel)
short recbuf[2*N+2*M]; // recording buffer
short recbufindex = 0;
short max_recbuf = 0;
short playback_scale = 1;
int state = 0;
short vclock_counter = 0; // virtual clock counter
short recbuf_start_clock = 0; // virtual clock counter for first sample in recording buffer
char r = 0;
short max_samp = 0;
float minz = 1E12;
float maxz = 0.0;

// highpass filter
short inindex = 0;  // input buffer index
float hpfout = 0.0;
short ii = 0;
short kk = 0;
//#define BL 21
//const float B[BL] = {
//    -0.0130241951,  0.02062470838,  0.02697116137, -0.01342090964, -0.04743985832,
//   -0.01095565408,  0.06797167659,  0.06939861178, -0.08312042803,  -0.3051353693,
//      0.588696897,  -0.3051353693, -0.08312042803,  0.06939861178,  0.06797167659,
//   -0.01095565408, -0.04743985832, -0.01342090964,  0.02697116137,  0.02062470838,
//    -0.0130241951
//};
#define BL 29
const float B[BL] = {
  -0.002215331886, 0.003098073881, 0.008297695778, -0.01649039052,-0.004188231193,
    0.01305069309,  0.02016889863, -0.01101524103, -0.03962349147, -0.00929670874,
     0.0632205382,   0.0657851845,  -0.0823623687,  -0.3033712804,   0.5898267031,
    -0.3033712804,  -0.0823623687,   0.0657851845,   0.0632205382, -0.00929670874,
   -0.03962349147, -0.01101524103,  0.02016889863,  0.01305069309,-0.004188231193,
   -0.01649039052, 0.008297695778, 0.003098073881,-0.002215331886
};

float inbuf[BL];

DSK6713_AIC23_CodecHandle hCodec;							// Codec handle
DSK6713_AIC23_Config config = DSK6713_AIC23_DEFAULTCONFIG;  // Codec configuration with default settings
// ------------------------------------------
// end of variables
// ------------------------------------------

double sin(double);
double cos(double);
interrupt void serialPortRcvISR(void);

void main()
{

	// set up the cosine and sin matched filters for searching
	// also initialize searching buffer
	for (i=0;i<M;i++){
		t = i*0.25;				// time
		y = cos(2*PI*t);		// cosine matched filter (double)
		mfc[i] = (float) y;		// cast and store
		y = sin(2*PI*t);		// sine matched filter (double)
		mfs[i] = (float) y;     // cast and store
		buf[i] = 0.0;             // clear searching buffer
		fbuf[i] = 0.0;
	}

	// initialize input buffer
	for (i=0;i<BL;i++)
		inbuf[i] = 0.0;

	// initialize clock buffer
	for (i=0;i<L;i++)
		clockbuf[i] = 0;

	// set up clock buffer to play modulated sinc centered at zero
	for (i=-N;i<=N;i++){
		x = i*BW;
		if (i!=0) {
			t = i*0.25;
			y = ((double) CLOCKAMPLITUDE)*cos(CARRIERKHZ*PI*t)*sin(PI*x)/(PI*x); // double
		}
		else {
			y = ((double) CLOCKAMPLITUDE);
		}
		j = i;
		if (j<0) {
			j += L; // wrap
		}
		clockbuf[j] = (short) y;
	}



	DSK6713_init();		// Initialize the board support library, must be called first
	DSK6713_LED_init(); // initialize LEDs
    hCodec = DSK6713_AIC23_openCodec(0, &config);	// open codec and get handle

	// Configure buffered serial ports for 32 bit operation
	// This allows transfer of both right and left channels in one read/write
	MCBSP_FSETS(SPCR1, RINTM, FRM);
	MCBSP_FSETS(SPCR1, XINTM, FRM);
	MCBSP_FSETS(RCR1, RWDLEN1, 32BIT);
	MCBSP_FSETS(XCR1, XWDLEN1, 32BIT);

	// set codec sampling frequency
	DSK6713_AIC23_setFreq(hCodec, DSK6713_AIC23_FREQ_16KHZ);

	// interrupt setup
	IRQ_globalDisable();			// Globally disables interrupts
	IRQ_nmiEnable();				// Enables the NMI interrupt
	IRQ_map(IRQ_EVT_RINT1,15);		// Maps an event to a physical interrupt
	IRQ_enable(IRQ_EVT_RINT1);		// Enables the event
	IRQ_globalEnable();				// Globally enables interrupts

	while(1)						// main loop
	{
	}
}

interrupt void serialPortRcvISR()
{
	union {Uint32 combo; short channel[2];} temp;

	temp.combo = MCBSP_read(DSK6713_AIC23_DATAHANDLE);
	// Note that right channel is in temp.channel[0]
	// Note that left channel is in temp.channel[1]

	// keep track of largest sample for diagnostics
	if (temp.channel[0]>max_samp)
		max_samp = temp.channel[0];

	// filter
	inbuf[inindex] = (float) temp.channel[0];  // right channel

	if (CLOCKOUTPUTCHANNEL==0) {
		// superimposed clocks and sync pulses, need to use HPF to avoid false detections
		hpfout = 0.0;
		ii = 0;
		for (kk=inindex;kk>=0;kk--){
			hpfout += inbuf[kk]*B[ii];
			ii++;
		}
		for (kk=BL-1;kk>inindex;kk--){
			hpfout += inbuf[kk]*B[ii];
			ii++;
		}
	}
	else {
		// clocks on left output, no need for HPF
		hpfout = (float) temp.channel[0];  // right channel
	}

	if (state==0) {  // SEARCHING STATE

        // put sample in searching buffer
		fbuf[bufindex] = hpfout;  // right channel
		buf[bufindex] = (float) temp.channel[0];  // unfiltered right channel
		if (PASSTHROUGH!=1) {
			temp.channel[0] = 0;  // can comment this out for debugging at home
		}

	    // compute incoherent correlation on filtered buffer
	    zc = 0;
	    zs = 0;
		for(i=0;i<M;i++) {
			zc += mfc[i]*fbuf[i];
			zs += mfs[i]*fbuf[i];
		}
		z = zc*zc+zs*zs;

		if (z>T1) {  				// threshold exceeded?
			if (z<minz)
			minz = z;				// diagnostic
			if (z>maxz)
				maxz = z;			// diagnostic
			max_recbuf = 0;
			state = 1; 				// enter "recording" state (takes effect in next interrupt)
		    DSK6713_LED_toggle(0);  // toggle LED for diagnostics
			recbufindex = M;		// start recording new samples at position M
			j = bufindex;			//
			for (i=0;i<M;i++){  	// copy samples from buf to first M elements of recbuf
				j++;   				// the first time through, this puts us at the oldest sample
				if (j>=M)
					j=0;
				//recbuf[i] = (short) buf[j];  // Oct 28 bad idea
				recbuf[i] = (short) fbuf[j];  // DRB: filtered buffer copied now to avoid echoing back clock pulses
				fbuf[j] = 0;  		// clear out searching buffer to avoid false trigger
				buf[j] = 0;  		// clear out searching buffer to avoid false trigger
			}
		}
		else {
			// increment and wrap pointer
		    bufindex++;
		    if (bufindex>=M)
		    	bufindex = 0;
		}

	}
	else if (state==1) { // RECORDING STATE

		// put sample in recording buffer
		//recbuf[recbufindex] = temp.channel[0];  // right channel
		recbuf[recbufindex] = (short) hpfout;  // filtered right channel
		if (PASSTHROUGH!=1) {
			temp.channel[0] = 0;  // can comment this out for debugging at home
		}
		if (abs(recbuf[recbufindex])>max_recbuf) {
			max_recbuf = abs(recbuf[recbufindex]); // keep track of largest sample for diagnostics
		}
		recbufindex++;
		if (recbufindex>=(2*N+2*M)) {
			state = 2;  		// buffer is full (stop recording and start counting samples to next clock tick)
			DSK6713_LED_toggle(0); // toggle LED for diagnostics
			recbufindex = 0; 	// shouldn't be necessary
			wait_count = 0;     // reset the wait counter
			if (max_recbuf<2000)
				playback_scale = 0;  // don't send response (signal was too weak)
			else if (max_recbuf<4000)
				playback_scale = 8;  // reply and scale by 8 (DRB now 4 to avoid overflow)
			else if (max_recbuf<8000)
				playback_scale = 4;  // reply and scale by 4 (DRB now 2 to avoid overflow)
			else if (max_recbuf<16000)
			playback_scale = 2;  // reply and scale by 2 (DRB now 1 tp avoid overflow)
			else
				playback_scale = 1;  // no scaling
		}
	}
	else if (state==2) { // WAITING STATE (counting up samples until clock tick)
		if (PASSTHROUGH!=1) {
			temp.channel[0] = 0;  // can comment this out for debugging at home
		}
		wait_count++;
		// only way out of this state is when the virtual clock rolls over
	}
	else if (state==3) { // WAITING STATE (counting down samples until reply)
		if (PASSTHROUGH!=1) {
			temp.channel[0] = 0;  // can comment this out for debugging at home
		}
		wait_count--;
		if (wait_count<HPFGROUPDELAY) {
			state = 4;  // done waiting, time to play back rec_buffer in reverse
			recbufindex = 2*N+2*M;
		}
	}
	else if (state==4) { // RESPONSE STATE (play back recording buffer in reverse)
		recbufindex--;
		if (recbufindex>=0) {
			temp.channel[0] = playback_scale*recbuf[recbufindex];  // response always on right channel
		}
		else
		{
			state = 0;  // go back to searching
		}
	}

	// temp.channel[0] right channel has been previously set
	temp.channel[1] = 0;	// clear left channel
	temp.channel[CLOCKOUTPUTCHANNEL] += clockbuf[vclock_counter];  // master clock signal (always played)
//	temp.channel[0] = temp.channel[1];  // XXX DIAGNOSTIC FOR OUTPUT OFFSET TEST
//	temp.channel[0] = clockbuf[vclock_counter];  // XXX DIAGNOSTIC


	MCBSP_write(DSK6713_AIC23_DATAHANDLE, temp.combo); // output L/R channels

	// update virtual clock (cycles from 0 to L-1)
	vclock_counter++;
	if (vclock_counter>=L) {
		vclock_counter = 0; // clock tick occurred, wrap
		if (state==2) {
			state = 3;
		}
	}

	inindex++;
	if (inindex>BL-1)
		inindex = 0;

}

