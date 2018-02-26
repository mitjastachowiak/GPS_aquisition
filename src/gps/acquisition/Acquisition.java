package gps.acquisition;

import java.nio.BufferOverflowException;

public class Acquisition {
  private int dopplerShift = Integer.MIN_VALUE;
  private int codeShift = Integer.MIN_VALUE;
  private int sampleCount = 0;
  private int codeCount = 0;
  
  // working vectors
  private ComplexVec Xin;   // samples
  private ComplexVec C;     // codes
  private ComplexVec Xfd;

  // algorithm parameters
  private final float fs = 400000; // hz
  private final float fstep = 1000; // hz
  private final float fmax = 5000; // hz
  private final float fmin = -5000; // hz
  private final float sn_threshold = 0.015f;
  public static final boolean debug = true;
  
  public Acquisition(int nrOfSamples){
    Xin = new ComplexVec(nrOfSamples);
    C = new ComplexVec(nrOfSamples);
    Xfd = new ComplexVec(nrOfSamples);
    ComplexVec.initF(nrOfSamples);
  }
  
  static void log (String str) {
    if (debug) System.out.println(str);
  }
  
  public void enterSample(float real, float imag){
    if (debug && sampleCount == Xin.vec.length) throw new BufferOverflowException();
    Xin.vec[sampleCount].r = real;
    Xin.vec[sampleCount].i = imag;
    sampleCount++;
  }
  
  public void enterCode(float real, float imag){
    if (debug && codeCount == C.vec.length) throw new BufferOverflowException();
    C.vec[codeCount].r = real;
    C.vec[codeCount].i = imag;
    codeCount++;
  }
  
  private void computeXfd (int N, float fd) {
    for (int n = 0; n < N; n++) {
      Xfd.vec[n].takeFrom(Xin.vec[n]);
      Xfd.vec[n].mulBy(Complex.exp(new Complex(0, -2*(float)Math.PI*fd*n/fs)));
    }
  }
  
  public boolean startAcquisition() {
    log("");
    log("== Start Acquisition ==");
    if (debug && (codeCount != sampleCount || codeCount != Xin.vec.length)) throw new IllegalArgumentException("Required ammount of samples/codes was not given!");
    int m = (int)((fmax-fmin)/fstep) + 1;
    int N = Xin.vec.length;
    log("fft(C):");
    C.fft();
    C.adj();
    // search maximum in R (compute R on the fly)
    log("search max:");
    Max max = new Max();
    for (int nf = 0; nf < m; nf++) {
      float fd = fmin + nf*fstep;
      log("  == fd="+fd+" ==");
      // compute Xfd
      log("  compute Xfd...");
      computeXfd(N, fd);
      // compute rfd*N
      log("  fft(Xfd):");
      Xfd.fft();
      log("  compute rfd...");
      ComplexVec rfd = Xfd;
      rfd.mul(C);
      log("  invfft(rfd):");
      rfd.invfft();
      // search for maximum
      log("  update max...");
      max.searchIn(rfd, nf);
    }
    log("compute doppler and codeshift...");
    dopplerShift = (max.nf - (int)(m/2)) * (int)fstep;
    codeShift = max.i;
    // compute gamma = signal to noise
    log("compute gamma...");
    double Pin = 0;
    for (int k = 0; k < N; k++) Pin += Xin.vec[k].abs2();
    double gamma = max.m / Pin; // the original formular would be max / (Pin/N), but rfd is *N here.
    return gamma > sn_threshold;
  }
  
  public int getDopplerverschiebung () { return dopplerShift; }
  
  public int getCodeVerschiebung() { return codeShift;  }




  private static class Max {
    public float m = 0;
    public int nf = 0;
    public int i = 0;
    
    private void searchCompare (float abs, int i, int nf) {
      if (abs > this.m) {
        this.m = abs;
        this.nf = nf;
        this.i = i;
      }
    }
    public void searchIn (ComplexVec v, int nf) {
      for (int i = 0; i < v.vec.length; i++) {
        float abs = v.vec[i].abs2();
        searchCompare(abs, i, nf);
      }
    }
  }




  private static class Complex {
    public float r;
    public float i;
    public Complex (float r, float i) {
      this.r = r;
      this.i = i;
    }
    public Complex (Complex c) {
      this.r = c.r;
      this.i = c.i;
    }
    
    public void takeFrom (Complex c) { this.r = c.r; this.i = c.i; }
      
    public Complex add (Complex c) { return new Complex(this.r+c.r, this.i+c.i); }
    public void incBy (Complex c) { this.r += c.r; this.i += c.i; }
    public Complex sub (Complex c) { return new Complex(this.r-c.r, this.i-c.i); }
      
    public Complex mul (Complex c) { return new Complex(this.r*c.r - this.i*c.i, this.r*c.i + this.i*c.r); }
    public void mulBy (Complex c) { float _r = this.r; this.r = this.r*c.r - this.i*c.i; this.i = _r*c.i + this.i*c.r; }
      
    public void divBy (float real) { this.r /= real; this.i /= real; }
      
    public void conj () { this.i = -this.i; }
      
    public float abs2 () { return this.r*this.r + this.i*this.i; }
      
    static public Complex exp (Complex c) {
      float e = (float)Math.exp(c.r);
      return new Complex(e*(float)Math.cos(c.i), e*(float)Math.sin(c.i));
    }
  }





  private static class ComplexVec {
    public final Complex[] vec;
    private static Complex[] f; // helper array for dft and fft-splitting, initialized outside (minimum length = nSamples/2 if nSamples even, nSamples otherwise)
    private static int nFixed;
    
    public static void initF(int n) {
      nFixed = n;
      if (n % 4 == 0) n = n/2; // fft-seperation needs at least n/2
      if (n % 2 != 0) n = 2*n; // dft on highest layer requires 2 times the array size
      f = new Complex[n];
      for (int i = 0; i < n; i++) f[i] = new Complex(0, 0);
    }
    
    public ComplexVec (int size) {
      vec = new Complex[size];
      for (int i = 0; i < size; i++) vec[i] = new Complex(0, 0);
    }
      
    public void mul (ComplexVec v) { // multiplies each element in this vector with the related element in v
      if (debug && this.vec.length != v.vec.length) throw new IllegalArgumentException("Vectors don't have equal size!");
      for (int k = 0; k < this.vec.length; k++) this.vec[k].mulBy(v.vec[k]);
    }
      
    public void div (float real) { // divides all elements in this vector by r
      for (int k = 0; k < this.vec.length; k++) this.vec[k].divBy(real);
    }
      
    public void adj () { // complex adjunks all elements in this vector
      for (int k = 0; k < this.vec.length; k++) this.vec[k].conj();
    }
      
    private void separate (int n_2, int offset) { // computes the fft-separation of a n-element part of this vector beginning at offset
      for(int i=0; i<n_2; i++) f[i].takeFrom(this.vec[i*2+1+offset]);
      for(int i=0; i<n_2; i++) this.vec[i+offset].takeFrom(this.vec[i*2+offset]);
      for(int i=0; i<n_2; i++) this.vec[i+n_2+offset].takeFrom(f[i]);
    }
    private void separateAll (int n_2) { // computes a fft-separation on all n-element parts of this vector
      for (int offset = 0; offset < this.vec.length; offset += n_2 + n_2) {
        log("    separate(n="+n_2+", offset="+offset+")");
        separate(n_2, offset);
      }
    }
    private void dftInner (int k, int n, int offset) { // inner loop of dft-function, just because double loops are not possible
      for (int i = 0; i < n; i++) {
        f[k].incBy(f[n + (i*k)%n].mul(this.vec[i + offset]));
      }
    }
    private void dft (int n, int offset) { // computes a simple dft a n-elements part of this vector beginning at offset
      for (int k = 0; k < n; k++) {
        f[k].r = 0;
        f[k].i = 0;
        dftInner(k, n, offset);
      }
      for (int k = 0; k < n; k++) this.vec[k + offset].takeFrom(f[k]);
    }
    private void merge (int n, int offset) { // computes the fft-merge of n-elements part of this vector beginning at offset
      for (int k = 0; k < n; k++) {
        Complex e = new Complex(this.vec[k + offset]);
        Complex o = new Complex(this.vec[k + n + offset]);
        Complex w = Complex.exp(new Complex(0, -1.0f * (float)Math.PI * k / n));
        this.vec[k + offset].takeFrom(e.add(w.mul(o)));
        this.vec[k + n + offset].takeFrom(e.sub(w.mul(o)));
      }
    }
    private void mergeAll (int n) { // computes fft-merges of all n-elements parts of this vector
      for (int offset = 0; offset < this.vec.length; offset += n + n) {
        log("    merge(n="+n+", offset="+offset+")");
        merge(n, offset);
      }
    }
    public void fft () {
      if (debug && this.vec.length != nFixed) throw new IllegalArgumentException("ComplexVec.initF was not called for this vector's size. Can't do fft!");
      // determine, how often n can be divided by 2
      int n = this.vec.length;
      int nMin = n;
      while (nMin % 2 == 0) nMin /= 2;
      // go town till nMin
      for (n = n/2; n >= nMin; n/=2) {
        separateAll(n);
      }
      // do dft for all unsplittable parts
      for (int i = 0; i < nMin; i++) {
        f[i+nMin].takeFrom(Complex.exp(new Complex(0, -2.0f*(float)Math.PI*i/nMin)));
      }
      if (nMin > 1) for (int offset = this.vec.length - nMin; offset >= 0; offset -= nMin) {
        log("    dft(n="+nMin+", offset="+offset+")");
        dft(nMin, offset);
      }
      // go up till this.vec.length
      for (n = nMin; n < this.vec.length; n *= 2) {
        mergeAll(n);
      }
    }
    
    public void invfft () {
      this.adj();
      this.fft();
      this.div(this.vec.length);
    }
    
    public void print () {
      for (int k = 0; k < this.vec.length; k++) System.out.println("   "+Double.toString(this.vec[k].r) + " + j" + Double.toString(this.vec[k].i));
    }
  }
}
