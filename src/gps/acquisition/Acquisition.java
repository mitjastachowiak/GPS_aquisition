package gps.acquisition;

import java.nio.BufferOverflowException;

import javax.naming.directory.InvalidAttributesException;

public class Acquisition {
	private int dopplerShift = Integer.MIN_VALUE;
	private int codeShift = Integer.MIN_VALUE;
	private int sampleCount = 0;
	private int codeCount = 0;
	private float[][] samples; // Xin
	private float[][] codes;   // C
	
	public Acquisition(int nrOfSamples){
		samples = new float[nrOfSamples][2];
		codes = new float[nrOfSamples][2];
	}
	
	public void enterSample(float real, float imag){
		if (sampleCount == samples.length) throw new BufferOverflowException();
		samples[sampleCount][0] = real;
		samples[sampleCount][1] = imag;
		sampleCount++;
	}
	
	public void enterCode(float real, float imag){
		if (codeCount == codes.length) throw new BufferOverflowException();
		codes[codeCount][0] = real;
		codes[codeCount][1] = imag;
		codeCount++;
	}
	
	public boolean startAcquisition() {
		if (sampleCount != samples.length || codeCount != samples.length) throw new IllegalArgumentException("Required ammount of samples/codes was not given!");

		return false;
	}
	
	public int getDopplerverschiebung() { return dopplerShift; }
	
	public int getCodeVerschiebung(){ return codeShift; }
	
	

}
