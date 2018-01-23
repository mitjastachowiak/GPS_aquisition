package gps.acquisition;

import java.nio.BufferOverflowException;

import javax.naming.directory.InvalidAttributesException;

public class Acquisition {
	private int dopplerShift = Integer.MIN_VALUE;
	private int codeShift = Integer.MIN_VALUE;
	private int sampleCount = 0;
	private int codeCount = 0;
	private ComplexVec Xin;   // samples
	private ComplexVec C;     // codes
	
	public Acquisition(int nrOfSamples){
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
        float fs;
        float fd;
        
        ComplexVec Xfd = new ComplexVec(Xin.vec.length);
        for (int n = 0; n < Xin.vec.length; n++) {
          Xfd.vec[n] = Xin.vec[n].mul(Complex.j().mul(fd*n*2*Math.PI/fs));
        }

		ComplexVec rfd = Xfd.dft().mul(C.adj().dft()).idft().div(Xfd.vec.length);
		
		return false;
	}
	
	public int getDopplerverschiebung() { return dopplerShift; }
	
	public int getCodeVerschiebung(){ return codeShift; }
	
	
  public static class Complex {
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
	
	static public Complex j () { return new Complex(0, 1); }
	static public Complex exp (Complex c) {
	  double e = Math.exp(c.r);
	  return new Complex(e*Math.cos(c.i), e*Math.sin(c.i));
    }
  }
  
  
  
  public static class ComplexVec {
	public Complex[] vec;
	public ComplexVec (int size) {
      vec = new Complex[size];
	}
	
	public ComplexVec mul (ComplexVec v) {
	  if (this.vec.length != v.vec.length) throw new IllegalArgumentException("Vactors don't have equal size!");
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
  }
}
