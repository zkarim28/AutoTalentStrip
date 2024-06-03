#ifndef __PITCHSHIFTER__
#define __PITCHSHIFTER__

//#include "IPlug_include_in_plug_hdr.h"
//#include "IPopupMenuControl.h"
#include "fftsetup.h"
#include "mayer_fft.c"
//#include "IKnobMultiControlText.h"
//#include "Scales.h"
#include <math.h>

#define L2SC (float)3.32192809488736218171

class PitchShifter //: public IPlug
{
public:
    
    //Added these two vars
    unsigned long originalSampleRate; //Replacement for GetSampleRate();
    unsigned long fs; // Sample rate
    
  PitchShifter()
  {
      //Previously, The Constructor was filled with IGraphics initalizations for each knob and
      init(fs);
  }

  ~PitchShifter()
  {
    fft_des(fmembvars);
    free(cbi);
    free(cbo);
    free(cbonorm);
    free(cbwindow);
    free(hannwindow);
    free(acwinv);
    free(frag);
    free(ffttime);
    free(fftfreqre);
    free(fftfreqim);
  };

  void Reset()
  {
      //TRACE;
//      IMutexLock lock(this);
      
//      unsigned long sr = GetSampleRate();
      unsigned long sr = originalSampleRate;
    
    if( fs != sr) init(sr);
  }

  void OnParamChange(int paramIdx)
  {
      //based on the graphics interface that we use
  }

    void ProcessDoubleReplacing(double** inputs, double** outputs, int nFrames)
    {
        double* in1 = inputs[0];
        double* out1 = outputs[0];
        double* out2 = outputs[1];
        float fPersist = glidepersist;

        aref = (float)440 * pow(2, fTune / 12);

        unsigned long N = cbsize;
        unsigned long Nf = corrsize;

        long int ti;
        long int ti2;
        long int ti3;
        float tf;
        float tf2;
        float tf3;

        for (int s = 0; s < nFrames; ++s, ++in1, ++out1, ++out2)
        {
            // load data into circular buffer
            tf = (float)*in1;
            cbi[cbiwr] = tf;
            cbiwr++;
            if (cbiwr >= N) {
                cbiwr = 0;
            }

            // ********************
            // * Low-rate section *
            // ********************

            // Every N/noverlap samples, run pitch estimation / correction code
            if ((cbiwr) % (N / noverlap) == 0) {

                // ---- Obtain autocovariance ----

                // Window and fill FFT buffer
                ti2 = (long)cbiwr;
                for (ti = 0; ti < (long)N; ti++) {
                    ffttime[ti] = (float)(cbi[(ti2 - ti) % N] * cbwindow[ti]);
                }

                // Calculate FFT
                fft_forward(fmembvars, ffttime, fftfreqre, fftfreqim);

                // Remove DC
                fftfreqre[0] = 0;
                fftfreqim[0] = 0;

                // Take magnitude squared
                for (ti = 1; ti < (long)Nf; ti++) {
                    fftfreqre[ti] = (fftfreqre[ti]) * (fftfreqre[ti]) + (fftfreqim[ti]) * (fftfreqim[ti]);
                    fftfreqim[ti] = 0;
                }

                // Calculate IFFT
                fft_inverse(fmembvars, fftfreqre, fftfreqim, ffttime);

                // Normalize
                for (ti = 1; ti < (long)N; ti++) {
                    ffttime[ti] = ffttime[ti] / ffttime[0];
                }
                ffttime[0] = 1;

                //  ---- END Obtain autocovariance ----


                //  ---- Calculate pitch and confidence ----

                // Calculate pitch period
                //   Pitch period is determined by the location of the max (biased)
                //     peak within a given range
                //   Confidence is determined by the corresponding unbiased height
                tf2 = 0;
                pperiod = pmin;
                for (ti = nmin; ti < (long)nmax; ti++) {
                    ti2 = ti - 1;
                    ti3 = ti + 1;
                    if (ti2 < 0) {
                        ti2 = 0;
                    }
                    if (ti3 > (long)Nf) {
                        ti3 = Nf;
                    }
                    tf = ffttime[ti];

                    if (tf > ffttime[ti2] && tf >= ffttime[ti3] && tf > tf2) {
                        tf2 = tf;
                        conf = tf * acwinv[ti];
                        pperiod = (float)ti / fs;
                    }
                }

                // Convert to semitones
                pitch = (float)-12 * log10((float)aref * pperiod) * L2SC;
                pitch = pitch;
                pperiod = pperiod;
                conf = conf;

                //  ---- END Calculate pitch and confidence ----


                //  ---- Determine pitch target ----

                // If voiced
                if (conf >= vthresh) {
                    // TODO: Scale sliders
                    // Determine pitch target
                    tf = -1;
                    tf2 = 0;
                    tf3 = 0;
                    for (ti = 0; ti < 12; ti++) {
                        switch (ti) {
                        case 0:
                            tf2 = fNotes[9];
                            break;
                        case 1:
                            tf2 = fNotes[10];
                            break;
                        case 2:
                            tf2 = fNotes[11];
                            break;
                        case 3:
                            tf2 = fNotes[0];
                            break;
                        case 4:
                            tf2 = fNotes[1];
                            break;
                        case 5:
                            tf2 = fNotes[2];
                            break;
                        case 6:
                            tf2 = fNotes[3];
                            break;
                        case 7:
                            tf2 = fNotes[4];
                            break;
                        case 8:
                            tf2 = fNotes[5];
                            break;
                        case 9:
                            tf2 = fNotes[6];
                            break;
                        case 10:
                            tf2 = fNotes[7];
                            break;
                        case 11:
                            tf2 = fNotes[8];
                            break;
                        }
                        tf2 = tf2 - (float)fabs((pitch - (float)ti) / 6 - 2 * floorf(((pitch - (float)ti) / 12 + 0.5)));
                        if (tf2 >= tf) {
                            tf3 = (float)ti;
                            tf = tf2;
                        }
                    }
                    ptarget = tf3;

                    // Glide persist
                    if (wasvoiced == 0) {
                        wasvoiced = 1;
                        tf = persistamt;
                        sptarget = (1 - tf) * ptarget + tf * sptarget;
                        persistamt = 1;
                    }

                    // Glide on circular scale
                    tf3 = (float)ptarget - sptarget;
                    tf3 = tf3 - (float)12 * floorf(tf3 / 12 + 0.5);
                    if (fGlide > 0) {
                        tf2 = (float)1 - pow((float)1 / 24, (float)N * 1000 / (noverlap * fs * fGlide));
                    }
                    else {
                        tf2 = 1;
                    }
                    sptarget = sptarget + tf3 * tf2;
                }
                // If not voiced
                else {
                    wasvoiced = 0;

                    // Keep track of persist amount
                    if (fPersist > 0) {
                        tf = pow((float)1 / 2, (float)N * 1000 / (noverlap * fs * fPersist));
                    }
                    else {
                        tf = 0;
                    }
                    persistamt = persistamt * tf; // Persist amount decays exponentially
                }
                // END If voiced

                //  ---- END Determine pitch target ----


                // ---- Determine correction to feed to the pitch shifter ----
                tf = sptarget - pitch; // Correction amount
                tf = tf - (float)12 * floorf(tf / 12 + 0.5); // Never do more than +- 6 semitones of correction
                if (conf < vthresh) {
                    tf = 0;
                }
                lrshift = fShift + fAmount * tf;  // Add in pitch shift slider


                // ---- Compute variables for pitch shifter that depend on pitch ---
                phincfact = (float)pow(2, lrshift / 12);
                if (conf >= vthresh) {  // Keep old period when unvoiced
                    phinc = (float)1 / (pperiod * fs);
                    phprd = pperiod * 2;
                }
            }
            // ************************
            // * END Low-Rate Section *
            // ************************


            // *****************
            // * Pitch Shifter *
            // *****************

            // TODO: Pre-filter with some kind of filter (maybe cheby2 or just svf)
            // TODO: Use cubic spline interpolation

            // IMPROVE QUALITY OF PITCH SHIFTER!
            // what is the glitch at "lAaAack"? probably pitch shifter

            //   Better snippet management
            //   Pre-filter
            //   Cubic spline interp
            // Pitch shifter (overlap-add, pitch synchronous)
            //   Note: pitch estimate is naturally N/2 samples old
            phasein = phasein + phinc;
            phaseout = phaseout + phinc * phincfact;

            //   If it happens that there are no snippets placed at the output, grab a new snippet!
            /*     if (cbonorm[((long int)cbord + (long int)(N/2*(1 - (float)1 / phincfact)))%N] < 0.2) { */
            /*       fprintf(stderr, "help!"); */
            /*       phasein = 1; */
            /*       phaseout = 1; */
            /*     } */

            //   When input phase resets, take a snippet from N/2 samples in the past
            if (phasein >= 1) {
                phasein = phasein - 1;
                ti2 = cbiwr - (long int)N / 2;
                for (ti = -((long int)N) / 2; ti < (long int)N / 2; ti++) {
                    frag[ti % N] = cbi[(ti + ti2) % N];
                }
            }

            //   When output phase resets, put a snippet N/2 samples in the future
            if (phaseout >= 1) {
                fragsize = fragsize * 2;
                if (fragsize >= N) {
                    fragsize = N;
                }
                phaseout = phaseout - 1;
                ti2 = cbord + N / 2;
                ti3 = (long int)(((float)fragsize) / phincfact);
                for (ti = -ti3 / 2; ti < (ti3 / 2); ti++) {
                    tf = hannwindow[(long int)N / 2 + ti * (long int)N / ti3];
                    cbo[(ti + ti2) % N] = cbo[(ti + ti2) % N] + frag[((int)(phincfact * ti)) % N] * tf;
                    cbonorm[(ti + ti2) % N] = cbonorm[(ti + ti2) % N] + tf;
                }
                fragsize = 0;
            }
            fragsize++;

            //   Get output signal from buffer
            tf = cbonorm[cbord];
            //   Normalize
            if (tf > 0.5) {
                tf = (float)1 / tf;
            }
            else {
                tf = 1;
            }
            tf = tf * cbo[cbord]; // read buffer
            tf = cbo[cbord];
            cbo[cbord] = 0; // erase for next cycle
            cbonorm[cbord] = 0;
            cbord++; // increment read pointer
            if (cbord >= N) {
                cbord = 0;
            }

            // *********************
            // * END Pitch Shifter *
            // *********************


            // Write audio to output of plugin
            // Mix (blend between original (delayed) =0 and shifted/corrected =1)
            *out1 = *out2 = (double)fMix * tf + (1 - fMix) * cbi[(cbiwr - N + 1) % N];
        }
    }

    
//  void ProcessDoubleReplacing(double** inputs, double** outputs, int nFrames)
//  {
//    double* in1 = inputs[0];
//    double* out1 = outputs[0];
//    double* out2 = outputs[1];
//    float fPersist = glidepersist;
//    
//    aref = (float)440*pow(2,fTune/12);
//          
//    unsigned long N = cbsize;
//    unsigned long Nf = corrsize;
//    
//    long int ti;
//    long int ti2;
//    long int ti3;
//    float tf;
//    float tf2;
//    float tf3;
//    
//    for (int s = 0; s < nFrames; ++s, ++in1, ++out1, ++out2)
//    {
//      tf = (float) *in1;
//      cbi[cbiwr] = tf;
//      cbiwr++;
//      if (cbiwr >= N) {
//        cbiwr = 0;
//      }
//
//      if ((cbiwr)%(N/noverlap) == 0) {
//        ti2 = (long) cbiwr;
//        for (ti=0; ti<(long)N; ti++) {
//          ffttime[ti] = (float)(cbi[(ti2-ti)%N]*cbwindow[ti]);
//        }
//        
//        fft_forward(fmembvars, ffttime, fftfreqre, fftfreqim);
//        fftfreqre[0] = 0;
//        fftfreqim[0] = 0;
//        
//        for (ti=1; ti< (long)Nf; ti++) {
//          fftfreqre[ti] = (fftfreqre[ti])*(fftfreqre[ti]) + (fftfreqim[ti])*(fftfreqim[ti]);
//          fftfreqim[ti] = 0;
//        }
//        
//        fft_inverse(fmembvars, fftfreqre, fftfreqim, ffttime);
//        for (ti=1; ti<(long)N; ti++) {
//          ffttime[ti] = ffttime[ti] / ffttime[0];
//        }
//        ffttime[0] = 1;
//
//        tf2 = 0;
//        pperiod = pmin;
//        for (ti=nmin; ti<(long)nmax; ti++) {
//          ti2 = ti-1;
//          ti3 = ti+1;
//          if (ti2<0) {
//            ti2 = 0;
//          }
//          if (ti3>(long)Nf) {
//            ti3 = Nf;
//          }
//          tf = ffttime[ti];
//          
//          if (tf>ffttime[ti2] && tf>=ffttime[ti3] && tf>tf2) {
//            tf2 = tf;
//            conf = tf*acwinv[ti];
//            pperiod = (float)ti/fs;
//          }
//        }
//        
//        pitch = (float) -12*log10((float)aref*pperiod)*L2SC;
//        pitch = pitch;
//        pperiod = pperiod;
//        conf = conf;
//        
//        if (conf>=vthresh) {
//          tf = -1;
//          tf2 = 0;
//          tf3 = 0;
//          for (ti=0; ti<12; ti++) {
//            switch (ti) {
//              case 0:
//                tf2 = fNotes[9];
//                break;
//              case 1:
//                tf2 = fNotes[10];
//                break;
//              case 2:
//                tf2 = fNotes[11];
//                break;
//              case 3:
//                tf2 = fNotes[0];
//                break;
//              case 4:
//                tf2 = fNotes[1];
//                break;
//              case 5:
//                tf2 = fNotes[2];
//                break;
//              case 6:
//                tf2 = fNotes[3];
//                break;
//              case 7:
//                tf2 = fNotes[4];
//                break;
//              case 8:
//                tf2 = fNotes[5];
//                break;
//              case 9:
//                tf2 = fNotes[6];
//                break;
//              case 10:
//                tf2 = fNotes[7];
//                break;
//              case 11:
//                tf2 = fNotes[8];
//                break;
//            }
//            tf2 = tf2 - (float)fabs( (pitch-(float)ti)/6 - 2*floorf(((pitch-(float)ti)/12 + 0.5)) );
//            if (tf2>=tf) {
//              tf3 = (float)ti;
//              tf = tf2;
//            }
//          }
//          ptarget = tf3;
//          
//          if (wasvoiced == 0) {
//            wasvoiced = 1;
//            tf = persistamt;
//            sptarget = (1-tf)*ptarget + tf*sptarget;
//            persistamt = 1;
//          }
//          
//          tf3 = (float)ptarget - sptarget;
//          tf3 = tf3 - (float)12*floorf(tf3/12 + 0.5);
//          if (fGlide>0) {
//            tf2 = (float)1-pow((float)1/24, (float)N * 1000/ (noverlap*fs*fGlide));
//          }
//          else {
//            tf2 = 1;
//          }
//          sptarget = sptarget + tf3*tf2;
//        }
//        else {
//          wasvoiced = 0;
//          if (fPersist>0) {
//            tf = pow((float)1/2, (float)N * 1000/ (noverlap*fs*fPersist));
//          }
//          else {
//            tf = 0;
//          }
//          persistamt = persistamt * tf;
//        }
//        
//        tf = sptarget - pitch;
//        tf = tf - (float)12*floorf(tf/12 + 0.5);
//        if (conf<vthresh) {
//          tf = 0;
//        }
//        lrshift = fShift + fAmount*tf;
//        
//        phincfact = (float)pow(2, lrshift/12);
//        if (conf>=vthresh) {
//          phinc = (float)1/(pperiod*fs);
//          phprd = pperiod*2;
//        }
//      }
//
//      phasein = phasein + phinc;
//      phaseout = phaseout + phinc*phincfact;
//      
//      if (phasein >= 1) {
//        phasein = phasein - 1;
//        ti2 = cbiwr - (long int)N/2;
//        for (ti=-((long int)N)/2; ti<(long int)N/2; ti++) {
//          frag[ti%N] = cbi[(ti + ti2)%N];
//        }
//      }
//
//      if (phaseout >= 1) {
//        fragsize = fragsize*2;
//        if (fragsize >= N) {
//          fragsize = N;
//        }
//        phaseout = phaseout - 1;
//        ti2 = cbord + N/2;
//        ti3 = (long int)(((float)fragsize) / phincfact);
//        for (ti=-ti3/2; ti<(ti3/2); ti++) {
//          tf = hannwindow[(long int)N/2 + ti*(long int)N/ti3];
//          cbo[(ti + ti2)%N] = cbo[(ti + ti2)%N] + frag[((int)(phincfact*ti))%N]*tf;
//          cbonorm[(ti + ti2)%N] = cbonorm[(ti + ti2)%N] + tf;
//        }
//        fragsize = 0;
//      }
//      fragsize++;
//      
//      tf = cbonorm[cbord];
//      if (tf>0.5) {
//        tf = (float)1/tf;
//      }
//      else {
//        tf = 1;
//      }
//      tf = tf*cbo[cbord];
//      tf = cbo[cbord];
//      cbo[cbord] = 0;
//      cbonorm[cbord] = 0;
//      cbord++;
//      if (cbord >= N) {
//        cbord = 0;
//      }
//      
//      *out1 = *out2 = (double) fMix*tf + (1-fMix)*cbi[(cbiwr - N + 1)%N];
//    }
//  }

  void init(unsigned long sr)
  {
      originalSampleRate = sr;
    unsigned long ti;
    
    fs = sr;
    aref = 440;
    
    if (fs >=88200) {
      cbsize = 4096;
    }
    else {
      cbsize = 2048;
    }
    corrsize = cbsize / 2 + 1;
    
    pmax = 1/(float)70;
    pmin = 1/(float)700;
    
    pperiod = pmax;
    
    nmax = (unsigned long)(fs * pmax);
    if (nmax > corrsize) {
      nmax = corrsize;
    }
    nmin = (unsigned long)(fs * pmin);
    
    cbi = (float*) calloc(cbsize, sizeof(float));
    cbo = (float*) calloc(cbsize, sizeof(float));
    cbonorm = (float*) calloc(cbsize, sizeof(float));
    
    cbiwr = 0;
    cbord = 0;
    
    hannwindow = (float*) calloc(cbsize, sizeof(float));
    for (ti=0; ti<cbsize; ti++) {
      hannwindow[ti] = -0.5*cos(2*M_PI*ti/(cbsize - 1)) + 0.5;
    }
    
    cbwindow = (float*) calloc(cbsize, sizeof(float));
    for (ti=0; ti<(cbsize / 2); ti++) {
      cbwindow[ti+cbsize/4] = -0.5*cos(4*M_PI*ti/(cbsize - 1)) + 0.5;
    }
    
    noverlap = 4;
    
    fmembvars = fft_con(cbsize);
    
    ffttime = (float*) calloc(cbsize, sizeof(float));
    fftfreqre = (float*) calloc(corrsize, sizeof(float));
    fftfreqim = (float*) calloc(corrsize, sizeof(float));
    
    acwinv = (float*) calloc(cbsize, sizeof(float));
    for (ti=0; ti<cbsize; ti++) {
      ffttime[ti] = cbwindow[ti];
    }
    fft_forward(fmembvars, cbwindow, fftfreqre, fftfreqim);
    for (ti=0; ti<corrsize; ti++) {
      fftfreqre[ti] = (fftfreqre[ti])*(fftfreqre[ti]) + (fftfreqim[ti])*(fftfreqim[ti]);
      fftfreqim[ti] = 0;
    }
    fft_inverse(fmembvars, fftfreqre, fftfreqim, ffttime);
    for (ti=1; ti<cbsize; ti++) {
      acwinv[ti] = ffttime[ti]/ffttime[0];
      if (acwinv[ti] > 0.000001) {
        acwinv[ti] = (float)1/acwinv[ti];
      }
      else {
        acwinv[ti] = 0;
      }
    }
    acwinv[0] = 1;
    
    lrshift = 0;
    ptarget = 0;
    sptarget = 0;
    wasvoiced = 0;
    persistamt = 0;
    
    glidepersist = 100;
    
    vthresh = 0.8;
    
    phprdd = 0.01;
    phprd = phprdd;
    phinc = (float)1/(phprd * fs);
    phincfact = 1;
    phasein = 0;
    phaseout = 0;
    frag = (float*) calloc(cbsize, sizeof(float));
    fragsize = 0;
  }

private:
  void SetScale();
  
  bool scales[8][12];
  
  // parameters
  float fMix;
  float fShift;
  float fTune;
  float fNotes[12];
  float fAmount;
  float fGlide;
  
  fft_vars* fmembvars; // member variables for fft routine
  
  unsigned long cbsize; // size of circular buffer
  unsigned long corrsize; // cbsize/2 + 1
  unsigned long cbiwr;
  unsigned long cbord;
  float* cbi; // circular input buffer
  float* cbo; // circular output buffer
  float* cbonorm; // circular output buffer used to normalize signal
  
  float* cbwindow; // hann of length N/2, zeros for the rest
  float* acwinv; // inverse of autocorrelation of window
  float* hannwindow; // length-N hann
  int noverlap;
  
  float* ffttime;
  float* fftfreqre;
  float* fftfreqim;
  
  // VARIABLES FOR LOW-RATE SECTION
  float aref; // A tuning reference (Hz)
  float pperiod; // Pitch period (seconds)
  float pitch; // Pitch (semitones)
  float pitchpers; // Pitch persist (semitones)
  float conf; // Confidence of pitch period estimate (between 0 and 1)
  float vthresh; // Voiced speech threshold
  
  float pmax; // Maximum allowable pitch period (seconds)
  float pmin; // Minimum allowable pitch period (seconds)
  unsigned long nmax; // Maximum period index for pitch prd est
  unsigned long nmin; // Minimum period index for pitch prd est
  
  float lrshift; // Shift prescribed by low-rate section
  int ptarget; // Pitch target, between 0 and 11
  float sptarget; // Smoothed pitch target
  int wasvoiced; // 1 if previous frame was voiced
  float persistamt; // Proportion of previous pitch considered during next voiced period
  float glidepersist;
  
  // VARIABLES FOR PITCH SHIFTER
  float phprd; // phase period
  float phprdd; // default (unvoiced) phase period
  float phinc; // input phase increment
  float phincfact; // factor determining output phase increment
  float phasein;
  float phaseout;
  float* frag; // windowed fragment of speech
  unsigned long fragsize; // size of fragment in samples
};

#endif
