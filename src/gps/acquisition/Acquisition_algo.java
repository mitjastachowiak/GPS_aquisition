package gps.acquisition;

import java.nio.BufferOverflowException;

import javax.naming.directory.InvalidAttributesException;

public class Acquisition_algo {
	private int dopplerShift = Integer.MIN_VALUE;
	private int codeShift = Integer.MIN_VALUE;
	private int sampleCount = 0;
	private int codeCount = 0;
	private ComplexVec Xin;   // samples
	private ComplexVec C;     // codes
	private final double fs = 400000; // hz
	private final double fstep = 1000; // hz
	private final double fmax = 5000; // hz
	private final double fmin = -5000; // hz
	private final double sn_threshold = 0.015;
	
	public Acquisition_algo(int nrOfSamples){
		Xin = new ComplexVec(nrOfSamples);
		C = new ComplexVec(nrOfSamples);
	}
	
	public void enterSample(float real, float imag){
		if (sampleCount == Xin.vec.length) throw new BufferOverflowException();
		Xin.vec[sampleCount] = new Complex(real, imag);
		sampleCount++;
	}
	
	public void enterCode(float real, float imag){
		if (codeCount == C.vec.length) throw new BufferOverflowException();
		C.vec[codeCount] = new Complex(real, imag);
		codeCount++;
	}
	
	public boolean startAcquisition() {
		if (sampleCount != Xin.vec.length || codeCount != Xin.vec.length) throw new IllegalArgumentException("Required ammount of samples/codes was not given!");
		
		// Compute R
		int m = (int)((fmax-fmin)/fstep) + 1;
        int N = Xin.vec.length;
        ComplexMat R = new ComplexMat(m, N);
        for (int nf = 0; nf < m; nf++) {
          double fd = fmin + nf*fstep;
          
          // compute Xfd
          ComplexVec Xfd = new ComplexVec(N);
          for (int n = 0; n < N; n++) {
            Xfd.vec[n] = Xin.vec[n].mul(Complex.exp(new Complex(0, -2*Math.PI*fd*n/fs)));
          }
          
          /*
          // Test: Xfd and the fourier-transforms are correct!
          System.out.println(" === fd = " + Double.toString(fd) + " === ");       
          System.out.println("Xfd = ");
          Xfd.print();
          System.out.println("Xfd_dft = ");
          Xfd.dft().print();                   // equal to 01_intermed_dft_X_fd.dat
          System.out.println("Xfd_idft = ");
          Xfd.dft().idft().print();            // equal to Xdf.print()
          */

          // compute rfd
		  ComplexVec rfd = Xfd.dft().mul(C.dft().adj()).idft().div(N);
		  R.setColumn(nf, rfd);
        }
        
        // find max
        ComplexMat.Max max = R.abs2Max();
        dopplerShift = (max.ix - (int)(m/2)) * (int)fstep;
        codeShift = max.iy;
        
        // compute gamma = signal to noise
        double Pin = 0;
        for (int k = 0; k < N; k++) Pin += Xin.vec[k].abs2();
        Pin /= N;
        double gamma = max.val / Pin;

		return gamma > sn_threshold;
	}
	
	public int getDopplerverschiebung() { return dopplerShift; }
	
	public int getCodeVerschiebung(){ return codeShift; }
	
	
  private static class Complex {
	public double r;
	public double i;
	public Complex (double r, double i) {
	  this.r = r;
	  this.i = i;
	}
	
	public Complex add (Complex c) { return new Complex(this.r+c.r, this.i*c.i); }
	public void incBy (Complex c) { this.r += c.r; this.i += c.i; }
	
	public Complex mul (Complex c) { return new Complex(this.r*c.r - this.i*c.i, this.r*c.i + this.i*c.r); }
	public Complex mul (double real) { return new Complex(this.r*real, this.i*real); }
	
	public Complex div (double real) { return new Complex(this.r/real, this.i/real); }
	
	public Complex adj () { return new Complex(this.r, -this.i); }
	
	public double abs2 () { return this.r*this.r + this.i*this.i; }
	
	static public Complex j () { return new Complex(0, 1); }
	static public Complex exp (Complex c) {
	  double e = Math.exp(c.r);
	  return new Complex(e*Math.cos(c.i), e*Math.sin(c.i));
    }
  }
  
  
  
  private static class ComplexVec {
	public Complex[] vec;
	public ComplexVec (int size) {
      vec = new Complex[size];
	}
	
	public ComplexVec mul (ComplexVec v) {
	  if (this.vec.length != v.vec.length) throw new IllegalArgumentException("Vectors don't have equal size!");
	  ComplexVec r = new ComplexVec(this.vec.length);
	  for (int k = 0; k < this.vec.length; k++) r.vec[k] = this.vec[k].mul(v.vec[k]);
	  return r;
	}
	
	public ComplexVec div (double real) {
	  ComplexVec r = new ComplexVec(this.vec.length);
	  for (int k = 0; k < this.vec.length; k++) r.vec[k] = this.vec[k].div(real);
	  return r;
	}
	
	public ComplexVec adj () {
	  ComplexVec a = new ComplexVec(this.vec.length);
	  for (int k = 0; k < this.vec.length; k++) a.vec[k] = this.vec[k].adj();
	  return a;
	}
	
	public ComplexVec dft () {
	  int N = this.vec.length;
	  ComplexVec f = new ComplexVec(N);
	  for (int k = 0; k < N; k++) {
	    f.vec[k] = new Complex(0, 0);
	    for (int i = 0; i < N; i++) {
	      f.vec[k].incBy(Complex.exp(new Complex(0, -2*Math.PI*i*k/N)).mul(this.vec[i]));
	    }
	  }
	  return f;
	}
	
	public ComplexVec idft () {
	  int N = this.vec.length;
	  ComplexVec a = new ComplexVec(N);
	  for (int k = 0; k < N; k++) {
	  	a.vec[k] = new Complex(0, 0);
	  	for (int i = 0; i < N; i++) {
	      a.vec[k].incBy(Complex.exp(new Complex(0, 2*Math.PI*i*k/N)).mul(this.vec[i]));
		}
	  	a.vec[k] = a.vec[k].div(N);
	  }
	  return a;
	}
	
	public void print () {
	  for (int k = 0; k < this.vec.length; k++) System.out.println("   "+Double.toString(this.vec[k].r) + " + j" + Double.toString(this.vec[k].i));
	}
  }
  
  
  
  private static class ComplexMat {
	public Complex[][] mat;
	public ComplexMat(int dimX, int dimY) {
      mat = new Complex[dimX][dimY];
	}
	
	public void setColumn(int column, ComplexVec v) {
	  if (column >= mat.length || v.vec.length != mat[column].length) throw new IllegalArgumentException("Vectorsize or column is not within matrix dimension!");
      for (int k = 0; k < v.vec.length; k++) mat[column][k] = v.vec[k];
	}
	
	public static class Max {
	  public double val;
	  public int ix;
	  public int iy;
	}
	public Max abs2Max () {
	  Max r = new Max();
	  for (int ix = 0; ix < mat.length; ix++) for (int iy = 0; iy < mat[ix].length; iy++) if (mat[ix][iy].abs2() > r.val) {
		r.val = mat[ix][iy].abs2();
		r.ix = ix;
		r.iy = iy;
	  }
	  return r;
	}
  }
}
