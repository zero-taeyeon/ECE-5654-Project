// rsx_analyze.c: S.W. Ellingson, Virginia Tech, 2023 Jun 09
// Reads raw data file obtained using the "rtl_sdr" program, and analyzes the data
// ---
// TO COMPILE: gcc -o rsx_analyze rsx_analyze.c -lm -lfftw3
// TO RUN:     ./rsx_analyze <filename> <n0>
//   <filename>: file to analyze
//   <n0>:       number of samples to skip (use to skip over samples at beginning of file)
// ---
// Min, max, mean for I,Q (sample values range -127..+128)                       --> stdout
// Min, max, mean dist in signal space between consecutive I,Q samples           --> stdout
// csv-formatted time domain (I,Q,mag,phase [deg])                               --> rsx_analyze_t.dat
// csv-formatted power spectrum, RSX_FFT_SIZE bins (bin,power[linear],power[dB]) --> rsx_analyze_f.dat
// csv-formatted time domain (I^2+Q^2 averaged over RSX_FFT_SIZE samples)        --> rsx_analyze_p.dat
// ---
// See end of this file for history.

#include <stdio.h>
//#include <stdlib.h>    /* for exit() */
#include <math.h>
#include <fftw3.h>

#define RSX_FFT_SIZE 2048
#define RSX_PI 3.1415927

int main ( int narg, char *argv[] ) {

  /*=================*/
  /*=== Variables ===*/
  /*=================*/

  unsigned char x[2*RSX_FFT_SIZE];
  float xs[RSX_FFT_SIZE];

  FILE *fp;   
  FILE *fpiq; 
  FILE *fpp; 
  
  /* FFTW variables */
  fftw_complex *fft_in, *fft_out;
  fftw_plan p;

  long int l;
  long int ll;

  int bFirst = 1;
  int ii, ii_previous;
  int qq, qq_previous;
  double s;

  long int nfft=0; /* number of FFT (sample blocks) processed */
  double w[RSX_FFT_SIZE]; /* window function */
 
  long int n0=0; /* first sample to be processed (use to skip samples at beginning of file) */
  int i_min, i_max;
  long int i_sum = 0.0;
  int q_min, q_max;
  long int q_sum = 0.0;
  double s_min, s_max;
  double s_sum = 0.0;

  long int nr;
  double rms_sum = 0;
  double i2q2 = 0.0; /* I^2+Q^2 */
  
  double i2q2a = 0.0; /* used to compute time average power */

  /* command line arguments */
  long int n=0;  /* number of samples processed */
  char filename[1024];

  /* Process command line arguments */
  sprintf(filename,"./capture.out");
  if (narg>1) {
    sscanf(argv[1],"%s",filename); 
    //} else {
    //  printf("FATAL: argument <dt> is required\n"); return 1; 
    }
  printf("rsx_analyze: filename = '%s'\n",filename);
  if (narg>2) {
    sscanf(argv[2],"%ld",&n0); 
    //} else {
    //  printf("FATAL: argument <dt> is required\n"); return 1; 
    }
  printf("rsx_analyze: n0 = %ld samples\n",n0);

  /* open file */
  fp = fopen(filename,"rb");
  if (!fp) {
    printf("rsx_analyze: FATAL: Unable to fopen() requested file\n");
    return 0;
    }

  /* set up FFTW */
  fft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * RSX_FFT_SIZE);
  fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * RSX_FFT_SIZE);
  p = fftw_plan_dft_1d( RSX_FFT_SIZE, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);

  /* initialize spectrogram */  
  for ( l=0; l<RSX_FFT_SIZE; l++ ) {
    xs[l] = 0.0;
    }

  /* setting up Blackman-Harris widow */
  /* https://en.wikipedia.org/wiki/Window_function */
  for ( l=0; l<RSX_FFT_SIZE; l++ ) {
    w[l] = 0.35875 - 0.48829*cos(2.0*RSX_PI*l/RSX_FFT_SIZE)
                   + 0.14128*cos(4.0*RSX_PI*l/RSX_FFT_SIZE)
                   - 0.01168*cos(6.0*RSX_PI*l/RSX_FFT_SIZE);
    }  

  /* initialize sample statistics */
  
  i_min = +129; i_max = -129; i_sum = 0;
  q_min = +129; q_max = -129; q_sum = 0;  
  s_min = 1.0e+10; s_max = -1.0e+10; s_sum = 0.0;

  /* open time-domain output files */
  fpiq = fopen("rsx_analyze_t.dat","w");
  fpp  = fopen("rsx_analyze_p.dat","w");
  
  /* prime the pump */
  fseek(fp,sizeof(x[0])*2*n0, SEEK_SET);   /* skip n0 samples at beginning of file */
  fread(x,sizeof(x[0]),2*RSX_FFT_SIZE,fp); /* get first batch */

  rms_sum = 0.0; /* this is for computing RMS over dataset */

  while (!feof(fp)) {

    nfft++;
    n += RSX_FFT_SIZE;

    i2q2a = 0.0;   /* this is for computing power over FFT length */

    for ( l=0; l<RSX_FFT_SIZE; l++ ) {
      ii = x[2*l  ] -127.0;
      qq = x[2*l+1] -127.0;
      
      /* process sample statistics */
      if (ii<i_min) { i_min = ii; }
      if (ii>i_max) { i_max = ii; }
      i_sum += ii;      
      if (qq<q_min) { q_min = qq; }
      if (qq>q_max) { q_max = qq; } 
      q_sum += qq;
      i2q2 = ((double)(ii*ii + qq*qq));
      rms_sum += i2q2;
      i2q2a   += i2q2;   
       
      if (!bFirst) {
	s = sqrt( (ii-ii_previous)*(ii-ii_previous) + (qq-qq_previous)*(qq-qq_previous) );
	//if (s>30.0) printf("%ld: %9.3e %+4d %+4d %+4d %+4d \n",n,s,ii,ii_previous,qq,qq_previous);		
        if (s<s_min) { s_min = s; }
        if (s>s_max) { s_max = s; }
        s_sum += s; 		   
	    } 
	  bFirst = 0;	 
      ii_previous = ii;
      qq_previous = qq; 
      
      fprintf(fpiq,"%+4d %+4d %7.3f %+8.3f\n",ii,qq,sqrt(ii*ii+qq*qq),atan2(qq,ii)*180.0/RSX_PI);
           
      } /* for l */		    

    fprintf(fpp,"%7.3f\n",10.0*log10(i2q2a/RSX_FFT_SIZE));

    /* loading into FFT input array while applying window */
    for ( l=0; l<RSX_FFT_SIZE; l++ ) {
      fft_in [l][0] = ( x[2*l  ] -127.0 ) * w[l];
      fft_in [l][1] = ( x[2*l+1] -127.0 ) * w[l]; 
      fft_out[l][0] = 0.0;  /* clearing output array, just in case */
      fft_out[l][1] = 0.0; 
      } /* for l */

    /* do that crazy FFT thing */	
    fftw_execute(p);

    /* copy result into spectrogram buffer, applying FFT shift as we go along */  
    for ( l=0; l<RSX_FFT_SIZE; l++ ) {
      ll = l + (RSX_FFT_SIZE/2); if (ll>=RSX_FFT_SIZE) { ll -= RSX_FFT_SIZE; } 
      xs[ll] += ( fft_out[l][0]*fft_out[l][0] + fft_out[l][1]*fft_out[l][1] )/(RSX_FFT_SIZE*RSX_FFT_SIZE);
      }

    /* read data */
    nr = fread(x,sizeof(x[0]),2*RSX_FFT_SIZE,fp);
    //printf("nr=%ld feof(fp)=%d\n",nr,feof(fp));

    } /* while (!feof(fp) */

  /* shut down FFTW */
  fftw_destroy_plan(p);
  fftw_free(fft_in); 
  fftw_free(fft_out);

  /* close input file */
  fclose(fp);

  /* close time-domain output files */
  fflush(fpiq); fclose(fpiq); printf("Time domain IQ data written to rsx_analyze_t.dat\n");
  fflush(fpp);  fclose(fpp);  printf("Time domain power data written to rsx_analyze_p.dat\n");
  
  /* report sample statistics */
  printf("%ld <-- Number of samples in file\n",n0+n+nr/2);
  printf("%ld <-- Number of samples processed\n",n); 
  printf("%+4d %+10.3e %+4d <-- I min, mean, max\n",i_min,((float)i_sum)/n,i_max);
  printf("%+4d %+10.3e %+4d <-- Q min, mean, max\n",q_min,((float)q_sum)/n,q_max);
  printf("%9.3e %9.3e %9.3e <-- dist. between adj. samples; min, mean, max\n",s_min,s_sum/(n-1),s_max);
  printf("%+10.3e <-- RMS sample magnitude\n",sqrt(rms_sum/n)); 

  /* write spectrogram, normalizing along way */ 
  fp = fopen("rsx_analyze_f.dat","w"); 
  for ( l=0; l<RSX_FFT_SIZE; l++ ) {
    //fprintf(fp,"%5ld %e %e\n",l,xs[l]/nfft,10.0*log10(xs[l]/nfft));
    fprintf(fp,"%+1.5f %e %e\n",((float)(l-(RSX_FFT_SIZE/2)))/RSX_FFT_SIZE,xs[l]/nfft,10.0*log10(xs[l]/nfft));
    }
  fflush(fp);
  fclose(fp);
  printf("Freq domain data written to rsx_analyze_f.dat\n");
  
  return 0;
  } /* main() */


//==================================================================================
//=== HISTORY ======================================================================
//==================================================================================

// rsx_analyze.c: S.W. Ellingson, Virginia Tech, 2023 Jun 09
// -- added time-averaged time-domain output
// rsx_analyze.c: S.W. Ellingson, Virginia Tech, 2019 Jul 24

//==================================================================================
//=== BELOW THIS LINE IS SCRATCH ===================================================
//==================================================================================



