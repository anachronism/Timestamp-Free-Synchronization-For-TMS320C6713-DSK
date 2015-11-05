/*************************************************************************
 *  Timestamp free synchronization
 *  Slave node code with coarse delay estimation
 *  DRB Apr 14, 2015
 *  Updated by Max Li Oct 12, 2015
 *  Updated by DRB October 24, 2015
 *    - Moved KF updates to ISR (and changed KF state units samples and samples/sample)
 *    - Fixed abs() bug in delay estimator (doesn't seem to make any difference)
 *  Updated again by DRB October 26, 2015
 *    - shifted clock pulse now superimposed on Line Out Right (at 1kHz carrier freq)
 *
 *  Plays out a modulated sinc pulse on the left channel every 4 clock periods

 *  States:
 *  State -1: warmup (needed since codec inputs seem to glitch on startup)
 *  State 0: searching
 *  State 1: recording
 *  State 2: coarse/fine delay estimation
 *
 *  TODO: Eliminate short offset. (DRB ?)
 *        Clean up unused constants and variables.
 *
 *************************************************************************/

// -------------------------------------
// DEFINES
// -------------------------------------
#define CHIP_6713 1

// uncomment one line below for corresponding M/S configuration (based on offsets from Google spreadsheet)
// units are in SAMPLES (not seconds)
//#define RECIPROCITY_CAL 5.92E-4  // A/B tests +37ns
//#define RECIPROCITY_CAL -9.60E-4  // B/A tests -60ns
//#define RECIPROCITY_CAL 4.32E-4  // A/C tests +27ns
//#define RECIPROCITY_CAL -7.76E-4  // C/A tests -48.5ns
#define RECIPROCITY_CAL -3.68E-4  // B/C tests -23ns
//#define RECIPROCITY_CAL 8.00E-6  // C/B tests 0.5ns
//#define RECIPROCITY_CAL 0 // no reciprocity calibration (default)

// clock pulse carrier kHz
// until October 2015, this was always 2
#define CARRIERKHZ 2

// which channel to output clock pulses
// until October 2015, this was always the left channel (0=right, 1=left)
// if this is set to 1, then the HPF is not run
#define CLOCKOUTPUTCHANNEL 1

// clock and sync pulse amplitude (should be no more than 16383 if pulses are added)
#define CLOCKAMPLITUDE 32000

// high pass filter *negative* group delay (HPF has group delay of 10) // -14
#define HPFGROUPDELAY 0

// length of searching window in samples
#define M 40

// threshold value for searching window
#define T1 100000

// sinc pulse normalized bandwidth
#define BW 0.0125
#define CBW 0.25 	//2kHz@8k Fs  carrier frequency

// 2*N+1 is the number of samples in the sinc function
#define N 560  // was 512

#define RECBUF_SIZE 1200 // was 1104 //2 * N + 2 * M

// virtual clock counter
#define L 4096
#define LOVER2 2048  // constant L/2

// modulated sinc pulse buffer length (mostly zeros)
#define LL 16384

// number of saves for testing/debugging
//#define SAVES 60
#define SAVES 200

// Kalman gains

#define K1 0.9 //Kalman gain coefficient 1 (used to be 0.95)
#define K2 0.00003 //Kalman gain coefficient 2 (used to be 0.45)
#define INHIBITSTATE2 10

#define W 1000000 //Number of entries in lookup tables
// define PI and INVPI
#define PI 3.14159265358979323846
#define INVPI 0.318309886183791

// -------------------------------------
// INCLUDES
// -------------------------------------
#include <stdio.h>
#include <c6x.h>
#include <csl.h>
#include <csl_mcbsp.h>
#include <csl_irq.h>
#include <math.h>

#include "dsk6713.h"
#include "dsk6713_aic23.h"
#include "DSK6713_LED.h"
#include "lookup_table.h"

// ------------------------------------------
// VARIABLES
// ------------------------------------------

// BUFFERS
float buf[M];       // search buffer
float fbuf[M];      // filtered search buffer
float mfc[M];		// in-phase correlation buffer
float mfs[M];       // quadrature correlation buffer
short sincpulsebuf[LL];  // sinc pulse buffer (left channel)
short clockbuf[L];  // clock buffer (right channel)
short clockbuf_shifted[L];  // double clock buffer (right channel)
float recbuf[2*N+2*M]; // recording buffer
float si[2*N+1];  // in-phase sinc pulse
float sq[2*N+1];  // quadrature sinc pulse

// COUNTERS, INDICES, FLAGS
short bufindex = 0;
short recbufindex = 0;
short save_index = 0;
int sincpulsebuf_counter = 0;  // sinc pulse buffer counter
short vclock_counter = 0; // virtual clock counter
short recbuf_start_clock = 0; // virtual clock counter for first sample in recording buffer
char hasSaved = 0; //Continues the circular nature of the buffers.
char allow_state2_update = 0;  // flag to stop KF from updating second state until first state is stable
char delay_est_done = 0;  // flag to tell ISR when delay estimates have been calculated
char new_observation = 0;  // flag to tell KF to use new observation

// Delay estimation variables
double clockoffset = 0.0;

// Used only in filling the shifted clock buffer the first time.
double fractionalShift = 0.0;
double t,x,y;

// Set up the state machine
short state = -1;  // start in state -1 to let codec settle

// Kalman filter and Just In Time buffer filling variables
double t_master_rightnow = 0.0;  // samples
double t_master_est = 0.0;  // master clock time estimates (samples)
double t_master_pred = 0.0;  // master clock time predictions (samples)
double r_master_est = 1.0;  // master clock rate estimates (samples)
double r_master_pred = 1.0;  // master clock rate predictions (samples)

// highpass filter
short inindex = 0;  // input buffer index
//const int BL = 21;
//const float B[21] = {
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

// Saves
short vclock_counter_save[SAVES];
int sincpulsebuf_counter_save[SAVES];
double clockoffset_save[SAVES];
double ytilde_save[SAVES];
double t_master_pred_save[SAVES];

// Codec
DSK6713_AIC23_CodecHandle hCodec;							// Codec handle
DSK6713_AIC23_Config config = DSK6713_AIC23_DEFAULTCONFIG;  // Codec configuration with default set_ISRings

// lookups table for cos and sinc
#pragma DATA_SECTION(cosbuf,".mydata")
far float cosbuf[2 * W];
float *sincbuf; //This is only a workaround of using the external memory.

// ------------------------------------------
// PROTOTYPES
// ------------------------------------------

double sin(double);
double cos(double);
interrupt void serialPortRcvISR(void);

void main()
{

	// local variables to make sure ISR doesn't mess with them
	int i = 0;
	int j = 0;
	int imax = 0;
	double z = 0.0;
	double zmax = 0.0;
	double zz = 0.0;
	//double zi = 0.0;
	double zq = 0.0;
	short cde = 0;			// coarse delay estimate (integer)
	double pcf = 0.0;		// phase correction factor
	double fde = 0.0;		// fine delay estimate

	// initialize the save buffers
	for (i=0;i<SAVES;i++) {
		ytilde_save[i] = 0.0;
		vclock_counter_save[i] = 0;
		sincpulsebuf_counter_save[i] = 0;
		t_master_pred_save[i] = 0.0;
		clockoffset_save[i] = 0.0;
	}

	// Initialize the sinc buffer's start address
	sincbuf = &(cosbuf[W]);

	// Initialize the trig tables
	for (i = 0; i < W; i++){
		cosbuf[i] = cos(2 * PI * i / W);
		zz = -(N + 1)*BW*PI + (2*(N + 1)*BW*PI)/W * i; // (zz should range from -N*BW*PI to N*BW*PI )
		if(zz != 0)
			sincbuf[i] = sin(zz) / zz;
		else
			sincbuf[i] = 1;
	}

	// set up the fractionally shifted buffers

	// set up the cosine and sin matched filters for searching
	// also initialize searching buffer
	for (i=0;i<M;i++){
		t = i*CBW;				// time
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

	// initialize clock buffers
	for (i=0;i<L;i++) {
		clockbuf[i] = 0;
		clockbuf_shifted[i] = 0;
	}

	// initialize sinc pulse buffer
	for (i=0;i<LL;i++)
		sincpulsebuf[i] = 0;

	// set up clock buffer and sinc pulse buffer
	// to play modulated sinc centered at zero
	for (i=-N;i<=N;i++){
		x = i*BW;
		if (i!=0) {
			t = i*CBW;
			y = ((double) CLOCKAMPLITUDE)*cos(2*PI*t)*sin(PI*	x)/(PI*x); // double
		}
		else {
			y = (double) CLOCKAMPLITUDE;
		}
		j = i;
		if (j<0) {
			j += L; // wrap
		}
		clockbuf[j] = (short) y;
		clockbuf_shifted[j] = (short)y;
		j = i;
		if (j<0) {
			j += LL; // wrap
		}
		sincpulsebuf[j] = (short) y;
	}

	// set up inphase and quadrature sinc pulses for coarse and fine delay estimators
	j = 0;
	for (i=-N;i<=N;i++){
		x = i*BW;
		if (i!=0) {
			t = i*CBW;
			si[j] = (float) (cos(2*PI*t)*sin(PI*x)/(PI*x));
			sq[j] = (float) (sin(2*PI*t)*sin(PI*x)/(PI*x));
		}
		else {
			si[j] = 1.0;
			sq[j] = 0.0;
		}
		j++;
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

	//////DSK6713_LED_toggle(3);	// toggle LED here for diagnostics

	while(1)						// main loop
	{

		if ((state==2)&&(delay_est_done==0)){  // TIME TO COMPUTE DELAY ESTIMATES

			DSK6713_LED_on(3);


			///DSK6713_LED_toggle(1);	// toggle LED here for diagnostics

			// compute coarse delay estimate (fixed DRB 20151024 - no longer uses abs(z))
			zmax = 0.0;				// maximum
			imax = 0;				// index of maximum

			for (i=0;i<(2*M-1);i++){  // lag index
				z = 0;
				for (j=0;j<(2*N+1);j++) {
					z += si[j]*recbuf[i+j];  // correlation at lag i
				}
				if (z>zmax) {
					zmax = z;  // store maximum
					imax = i;  // store index of maximum
				}
			}

			// catch error
			if (imax<0){
				printf("err1\n");
			    while (1) {}
			}

			// catch error
			if (imax>=2*M){
				printf("err2\n");
				while (1) {}
			}

			cde = recbuf_start_clock + imax + N; // coarse delay estimate (DRB: +N here because si is already shifted by N)
			// cde  is the number of samples elapsed since we launched the S->M sinc pulse

		    // compute fine delay estimate
			// no need to compute zi because zi will equal zmax
			//zi = 0.0;  // in phase
			zq = 0.0;  // quadrature
			for (j=0;j<(2*N+1);j++) {
				//zi += si[j]*recbuf[imax+j];  // correlation at lag imax
				zq += sq[j]*recbuf[imax+j];  // correlation at lag imax
			}

			// catch error
			if isnan(zmax){
				printf("err3\n");
				while (1) {}
			}

			// catch error
			if isnan(zq){
				printf("err4\n");
				while (1) {}
			}

			pcf = atan2(zq,zmax)*(2*INVPI); // [units of samples] assumes wc = pi/2
			fde = cde + pcf + ((double) (HPFGROUPDELAY));  // [units of samples]

			// compute actual clock offset
			// the value computed here is always non-negative and
			// represents the number of samples the master clock is ahead of the slave clock
			// (or the number of samples the slave clock is behind the master clock)
			// S->M sinc pulse was launched at time 0
			// M->S sinc pulse was received at time fde (should be positive)
			// to synchronize, we want slave clock ticks to appear at fde/2 + k*L for k=0,1,....
			clockoffset = fde*0.5;  // [units of samples]
			//while (clockoffset>((double) L))
			//	clockoffset = clockoffset - (double) L;  // DRB I don't think we want to do this

			// tell the ISR the calculations are done
			delay_est_done = 1;
			new_observation = 1;
			///DSK6713_LED_toggle(1);	// toggle LED here for diagnostics

			DSK6713_LED_off(3);

			}

	}  // while(1)
}  // void main

interrupt void serialPortRcvISR()
{
	union {Uint32 combo; short channel[2];} temp;
	short ii = 0;
	short jj = 0;
	double x_ISR, t_ISR;
	double ytilde = 0.0;	// innovation
	float hpfout = 0.0;
    float zz, zc, zs;

	temp.combo = MCBSP_read(DSK6713_AIC23_DATAHANDLE);
	// Note that right channel is in temp.channel[0]
	// Note that left channel is in temp.channel[1]

	// filter
	inbuf[inindex] = (float) temp.channel[0];  // right channel

	if (CLOCKOUTPUTCHANNEL==0) {
		// CLOCKOUTPUTCHANNEL = 0 means we are superimposing clocks on the right channel
		// and that we need to run the HPF to avoid false detections of sinc pulses
		hpfout = 0.0;  // not necessary but leave in
		ii = 0;
		for (jj=inindex;jj>=0;jj--){
			hpfout += inbuf[jj]*B[ii];
			ii++;
		}
		for (jj=BL-1;jj>inindex;jj--){
			hpfout += inbuf[jj]*B[ii];
			ii++;
		}
	} else
	{
		// CLOCKOUTPUTCHANNEL = 1 means the clock is going out on the left channel
		// the HPF is bypassed if we output clocks on the left channel
		hpfout = (float) temp.channel[0];
	}

	if (state==0) {  // SEARCHING STATE

		DSK6713_LED_on(0);  // state 0 = LED 0 on

        // put sample in searching buffer
		fbuf[bufindex] = hpfout;  // filtered right channel
		buf[bufindex] = (float) temp.channel[0];  // unfiltered right channel

	    // compute incoherent correlation
	    zc = 0;
	    zs = 0;
		for (ii=0;ii<M;ii++) {
			zc += mfc[ii]*fbuf[ii];
			zs += mfs[ii]*fbuf[ii];
		}
		zz = zc*zc+zs*zs;

		if (zz>T1) {  				// threshold exceeded?
			////t_estimation++;
			state = 1; 				// enter "recording" state (takes effect in next interrupt)
			DSK6713_LED_off(0); // turn off LED
			///DSK6713_LED_toggle(0);	// toggle LED here for diagnostics
			// record time of first sample (DRB fixed 4/19/2015: added +1)
			recbuf_start_clock = sincpulsebuf_counter-M+1; // should not be negative since we started counting when we launched the S->M sinc pulse
			recbufindex = M;		// start recording new samples at position M
			jj = bufindex;			//
			for (ii=0;ii<M;ii++){  	// copy samples from buf to first M elements of recbuf
				jj++;   				// the first time through, this puts us at the oldest sample
				if (jj>=M)
					jj=0;
				recbuf[ii] = (short) fbuf[jj];  // DRB: unfiltered buffer now
				fbuf[jj] = 0;       // clear out
				buf[jj] = 0;  		// clear out searching buffer to avoid false trigger
			}
		}
		else {
			// increment and wrap pointer
		    bufindex++;
		    if (bufindex>=M)
		    	bufindex = 0;
		}

	}  // end of searching state
	else if (state==1) { // RECORDING STATE

		DSK6713_LED_on(1);  // state 1 = LED 1 on

		// put sample in recording buffer
		//recbuf[recbufindex] = (float) temp.channel[0];  // unfiltered right channel
		recbuf[recbufindex] = (short) hpfout;  // filtered
		recbufindex++;

		if (recbufindex>=(RECBUF_SIZE)) { //
			state = 2;  		// buffer is full
			DSK6713_LED_off(1);  // LED 1 off
			delay_est_done = 0; // clear flag
			///DSK6713_LED_toggle(0);	// toggle LED here for diagnostics
			recbufindex = 0; 	// shouldn't be necessary
		}

	}  // end of recording state
	else if (state==2) { // CALCULATING DELAY ESTIMATES STATE
		DSK6713_LED_on(2);  // state 2 = LED 2 on

		if (delay_est_done==1) {  // are the delay estimates done calculating?

			state = 0; // next state
			DSK6713_LED_off(2);  // LED off
		}
		//state -= 2 * delay_est_done;
	}

	/**************** UPDATE THE KF EVERY ISR  *****************/

    if ((new_observation==0)||(t_master_est<=N+3)||(t_master_est>=(L-(N + 3)))) {

    	t_master_est = t_master_pred;  // samples
    	r_master_est = r_master_pred;  // samples/sample
    	t_master_pred = t_master_est + r_master_est; // samples
    	// r_master_pred doesn't change

    }
    else { // new observation of clockoffset and we are not in the middle of playing out a pulse

    	// innovation
    	t_master_rightnow = (((double) sincpulsebuf_counter) - clockoffset)*r_master_est;
    	ytilde =  t_master_rightnow - t_master_pred;  // samples

    	// need to wrap ytilde to -L/2 to L/2
    	while (ytilde>=LOVER2)
    		ytilde -= L;
    	while (ytilde<-LOVER2)
    		ytilde +=L;

    	// update KF estimates
    	t_master_est = t_master_pred + ((double) K1)*ytilde;  // filtered clock offset estimate
    	if (allow_state2_update==1)
	        r_master_est = r_master_pred + ((double) K2)*ytilde;  // filtered frequency offset estimate
    	else
    		r_master_est = r_master_pred;

	    // update KF predictions
    	t_master_pred = t_master_est + r_master_est;
    	r_master_pred = r_master_est;

    	// update saves for diagnostics
    	vclock_counter_save[save_index] = vclock_counter;
    	sincpulsebuf_counter_save[save_index] = sincpulsebuf_counter;
    	clockoffset_save[save_index] = clockoffset;
    	t_master_pred_save[save_index] = t_master_pred;
    	ytilde_save[save_index] = ytilde;
    	save_index++;
    	if (save_index==INHIBITSTATE2)
    		allow_state2_update = 1;  // ytilde is small enough now that we can update the frequency estimates
    	if (save_index>=SAVES)
    		save_index = 0;

    	// clear flag
    	new_observation = 0;

    }

    // wrap the master time estimate to [0,L) so we can generate
    // a modulated sinc pulse on master time
    while (t_master_est>=L) {
    	t_master_est -= L;
    	t_master_pred -= L;
    }


	/************JUST IN TIME BUFFER CALCULATION*****************/


	//Calculate sinc value.
	if(t_master_est <= N + 1 ){
		x_ISR = (t_master_est+RECIPROCITY_CAL)*BW;
		if (x_ISR==0.0) {
			clockbuf_shifted[vclock_counter] = (double) CLOCKAMPLITUDE;
		}
		else
		{
			t_ISR = (t_master_est+RECIPROCITY_CAL)*CBW;
			clockbuf_shifted[vclock_counter] = ((double) CLOCKAMPLITUDE) * cos_lut(CARRIERKHZ * PI * t_ISR) * sinc_lut(x_ISR);  // carrier freq now 1kHz
		}
	}
	else if (t_master_est >= L-(N + 1)){
			x_ISR = (t_master_est+RECIPROCITY_CAL-L)*BW;
			if (x_ISR==0.0) {
				clockbuf_shifted[vclock_counter] = (double) CLOCKAMPLITUDE;
			}
			else
			{
				t_ISR = (t_master_est+RECIPROCITY_CAL-L)*CBW;
				clockbuf_shifted[vclock_counter] = ((double) CLOCKAMPLITUDE) * cos_lut(CARRIERKHZ * PI * t_ISR) * sinc_lut(x_ISR); // carrier freq now 1kHz
			}
		}
	else
		clockbuf_shifted[vclock_counter] = 0;


	/*******************OUTPUT BUFFER****************************/
	if (state==-1) { // WARMUP STATE (NO OUTPUT)
		if (sincpulsebuf_counter>LL/2) {
			state = 0;
		}
		temp.channel[1] = 0;
		temp.channel[0] = 0;
	}
	else {
		//temp.channel[1] = clockbuf[vclock_counter];  // slave *unshifted* clock signal (for debug)
		temp.channel[0] = sincpulsebuf[sincpulsebuf_counter];  // this initiates the sinc pulse exchange with the master
		temp.channel[1] = 0;  // clear left channel
		temp.channel[CLOCKOUTPUTCHANNEL] += clockbuf_shifted[vclock_counter];  // superimpose the shifted clock buffer on the right channel (NEW)
	}

	MCBSP_write(DSK6713_AIC23_DATAHANDLE, temp.combo); // output L/R channels

	/*******************UPDATE TIMERS****************************/
	vclock_counter++;
	if (vclock_counter>=L) {
		vclock_counter = 0; // clock tick occurred, wrap
	}

	// update sinc pulse counter (cycles from 0 to LL-1)
	// this is for sinc pulses from the slave to the master
	// sinc pulses from the slave to the master have their peak at
	// sincpulsebuf_counter = 0 (this makes clock offset calculations simple)
	sincpulsebuf_counter++;
	if (sincpulsebuf_counter>=LL) {
		sincpulsebuf_counter = 0; // wrap
	}

	inindex++;
	if (inindex>BL-1)
		inindex = 0;

}
