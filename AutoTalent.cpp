#include "AutoTalent.h"
//#include "IPlug_include_in_plug_src.h"
//#include "IControl.h"
//#include "resource.h"


//functions that definitely need to be in our new library

AutoTalent::AutoTalent()
{
    //Previously, The Constructor was filled with IGraphics initalizations
}

AutoTalent::~AutoTalent(){
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

void AutoTalent::ProcessDoubleReplacing(double** inputs, double** outputs, int nFrames)
{
  // Mutex is already locked for us.
  
  double* in1 = inputs[0];
  double* out1 = outputs[0];
  double* out2 = outputs[1];
  // copy struct variables to local
  /*
   float fMix = fMix;
   float fShift = fShift;
   float fTune = fTune;
   float fA = fA;
   float fBb = fBb;
   float fB = fB;
   float fC = fC;
   float fDb = fDb;
   float fD = fD;
   float fEb = fEb;
   float fE = fE;
   float fF = fF;
   float fGb = fGb;
   float fG = fG;
   float fAb = fAb;
   float fGlide = fGlide;
   float fAmount = fAmount;
   */
  float fPersist = glidepersist;
  
  aref = (float)440*pow(2,fTune/12);
        
  unsigned long N = cbsize;
  unsigned long Nf = corrsize;
  //unsigned long fs = fs;
  /*
   float pmax = pmax;
   float pmin = pmin;
   unsigned long nmax = nmax;
   unsigned long nmin = nmin;
   
   float pperiod = pperiod;
   float pitch = pitch;
   float conf = conf;
   float aref = aref;
   */
  //
  
  long int ti;
  long int ti2;
  long int ti3;
  float tf;
  float tf2;
  float tf3;
  
  
  //double samplesPerBeat = GetSamplesPerBeat();
  //double samplePos = (double) GetSamplePos();
  
  for (int s = 0; s < nFrames; ++s, ++in1, ++out1, ++out2)
  {
    
    // load data into circular buffer
    tf = (float) *in1;
    cbi[cbiwr] = tf;
    cbiwr++;
    if (cbiwr >= N) {
      cbiwr = 0;
    }
    
    
    // ********************
    // * Low-rate section *
    // ********************
    
    // Every N/noverlap samples, run pitch estimation / correction code
    if ((cbiwr)%(N/noverlap) == 0) {
      
      
      // ---- Obtain autocovariance ----
      
      // Window and fill FFT buffer
      ti2 = (long) cbiwr;
      for (ti=0; ti<(long)N; ti++) {
        ffttime[ti] = (float)(cbi[(ti2-ti)%N]*cbwindow[ti]);
      }
      
      // Calculate FFT
      fft_forward(fmembvars, ffttime, fftfreqre, fftfreqim);
      
      // Remove DC
      fftfreqre[0] = 0;
      fftfreqim[0] = 0;
      
      // Take magnitude squared
      for (ti=1; ti< (long) Nf; ti++) {
        fftfreqre[ti] = (fftfreqre[ti])*(fftfreqre[ti]) + (fftfreqim[ti])*(fftfreqim[ti]);
        fftfreqim[ti] = 0;
      }
      
      // Calculate IFFT
      fft_inverse(fmembvars, fftfreqre, fftfreqim, ffttime);
      
      // Normalize
      for (ti=1; ti<(long)N; ti++) {
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
      for (ti=nmin; ti<(long)nmax; ti++) {
        ti2 = ti-1;
        ti3 = ti+1;
        if (ti2<0) {
          ti2 = 0;
        }
        if (ti3>(long)Nf) {
          ti3 = Nf;
        }
        tf = ffttime[ti];
        
        if (tf>ffttime[ti2] && tf>=ffttime[ti3] && tf>tf2) {
          tf2 = tf;
          conf = tf*acwinv[ti];
          pperiod = (float)ti/fs;
        }
      }
      
      // Convert to semitones
      pitch = (float) -12*log10((float)aref*pperiod)*L2SC;
      pitch = pitch;
      pperiod = pperiod;
      conf = conf;
      
      //  ---- END Calculate pitch and confidence ----
      
      
      //  ---- Determine pitch target ----
      
      // If voiced
      if (conf>=vthresh) {
        // TODO: Scale sliders
        // Determine pitch target
        tf = -1;
        tf2 = 0;
        tf3 = 0;
        for (ti=0; ti<12; ti++) {
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
          /*       if (ti==ptarget) { */
          /*         tf2 = tf2 + 0.01; // add a little hysteresis */
          /*       } */
          tf2 = tf2 - (float)fabs( (pitch-(float)ti)/6 - 2*floorf(((pitch-(float)ti)/12 + 0.5)) ); // like a likelihood function
          if (tf2>=tf) {                                                                           // that we're maximizing
            tf3 = (float)ti;                                                                       // to find the target pitch
            tf = tf2;
          }
        }
        ptarget = tf3;
        
        // Glide persist
        if (wasvoiced == 0) {
          wasvoiced = 1;
          tf = persistamt;
          sptarget = (1-tf)*ptarget + tf*sptarget;
          persistamt = 1;
        }
        
        // Glide on circular scale
        tf3 = (float)ptarget - sptarget;
        tf3 = tf3 - (float)12*floorf(tf3/12 + 0.5);
        if (fGlide>0) {
          tf2 = (float)1-pow((float)1/24, (float)N * 1000/ (noverlap*fs*fGlide));
        }
        else {
          tf2 = 1;
        }
        sptarget = sptarget + tf3*tf2;
      }
      // If not voiced
      else {
        wasvoiced = 0;
        
        // Keep track of persist amount
        if (fPersist>0) {
          tf = pow((float)1/2, (float)N * 1000/ (noverlap*fs*fPersist));
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
      tf = tf - (float)12*floorf(tf/12 + 0.5); // Never do more than +- 6 semitones of correction
      if (conf<vthresh) {
        tf = 0;
      }
      lrshift = fShift + fAmount*tf;  // Add in pitch shift slider
      
      
      // ---- Compute variables for pitch shifter that depend on pitch ---
      phincfact = (float)pow(2, lrshift/12);
      if (conf>=vthresh) {  // Keep old period when unvoiced
        phinc = (float)1/(pperiod*fs);
        phprd = pperiod*2;
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
    phaseout = phaseout + phinc*phincfact;
    
    //   If it happens that there are no snippets placed at the output, grab a new snippet!
    /*     if (cbonorm[((long int)cbord + (long int)(N/2*(1 - (float)1 / phincfact)))%N] < 0.2) { */
    /*       fprintf(stderr, "help!"); */
    /*       phasein = 1; */
    /*       phaseout = 1; */
    /*     } */
    
    //   When input phase resets, take a snippet from N/2 samples in the past
    if (phasein >= 1) {
      phasein = phasein - 1;
      ti2 = cbiwr - (long int)N/2;
      for (ti=-((long int)N)/2; ti<(long int)N/2; ti++) {
        frag[ti%N] = cbi[(ti + ti2)%N];
      }
    }
    
    //   When output phase resets, put a snippet N/2 samples in the future
    if (phaseout >= 1) {
      fragsize = fragsize*2;
      if (fragsize >= N) {
        fragsize = N;
      }
      phaseout = phaseout - 1;
      ti2 = cbord + N/2;
      ti3 = (long int)(((float)fragsize) / phincfact);
      for (ti=-ti3/2; ti<(ti3/2); ti++) {
        tf = hannwindow[(long int)N/2 + ti*(long int)N/ti3];
        cbo[(ti + ti2)%N] = cbo[(ti + ti2)%N] + frag[((int)(phincfact*ti))%N]*tf;
        cbonorm[(ti + ti2)%N] = cbonorm[(ti + ti2)%N] + tf;
      }
      fragsize = 0;
    }
    fragsize++;
    
    //   Get output signal from buffer
    tf = cbonorm[cbord];
    //   Normalize
    if (tf>0.5) {
      tf = (float)1/tf;
    }
    else {
      tf = 1;
    }
    tf = tf*cbo[cbord]; // read buffer
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
    *out1 = *out2 = (double) fMix*tf + (1-fMix)*cbi[(cbiwr - N + 1)%N];
  }
  
}



void AutoTalent::init(unsigned long sr)
{
  
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
  
  pmax = 1/(float)70;  // max and min periods (ms)
  pmin = 1/(float)700; // eventually may want to bring these out as sliders
  
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
  
  // Standard raised cosine window, max height at N/2
  hannwindow = (float*) calloc(cbsize, sizeof(float));
  for (ti=0; ti<cbsize; ti++) {
    hannwindow[ti] = -0.5*cos(2*PI*ti/(cbsize - 1)) + 0.5;
  }
  
  // Generate a window with a single raised cosine from N/4 to 3N/4
  cbwindow = (float*) calloc(cbsize, sizeof(float));
  for (ti=0; ti<(cbsize / 2); ti++) {
    cbwindow[ti+cbsize/4] = -0.5*cos(4*PI*ti/(cbsize - 1)) + 0.5;
  }
  
  noverlap = 4;
  
  fmembvars = fft_con(cbsize);
  
  ffttime = (float*) calloc(cbsize, sizeof(float));
  fftfreqre = (float*) calloc(corrsize, sizeof(float));
  fftfreqim = (float*) calloc(corrsize, sizeof(float));
  
  
  // ---- Calculate autocorrelation of window ----
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
  // ---- END Calculate autocorrelation of window ----
  
  lrshift = 0;
  ptarget = 0;
  sptarget = 0;
  wasvoiced = 0;
  persistamt = 0;
  
  glidepersist = 100; // 100 ms glide persist
  
  vthresh = 0.8;  //  The voiced confidence (unbiased peak) threshold level
  
  // Pitch shifter initialization
  phprdd = 0.01; // Default period
  phprd = phprdd;
  phinc = (float)1/(phprd * fs);
  phincfact = 1;
  phasein = 0;
  phaseout = 0;
  frag = (float*) calloc(cbsize, sizeof(float));
  fragsize = 0;
}

void AutoTalent::Reset()
{
  TRACE;
  IMutexLock lock(this);
  
  unsigned long sr = GetSampleRate();
  
  if( fs != sr) init(sr);
}

void AutoTalent::OnParamChange(int paramIdx)
{
    //big case statement to implement using JUCE
    //Get the param and adjust accordinly
    
    //***divided every single parameter by 100.0 except for
    //kShift, kRoot, and kScale***
}

void AutoTalent::SetScale(){
//  int sc[12];
//  mScales.makeScale(GetParam(kRoot)->Value(), GetParam(kScale)->Value(), sc);
//  
//  for (int i = 0; i< 12; i++) {
//    Keys[i]->SetValueFromPlug(sc[i]);
//    Keys[i]->SetDirty(true);
//  }
    
    //in the case that we use scales to easily pitch bend to
    //Here was the previous implementation
}
