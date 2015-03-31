package io.compgen.cgseq.variant;

public class VariantResults {
	public final String majorCall;
	public final String minorCall;
	
	public final int majorCount;
	public final int minorCount;
	
	public final int bgCount;
	
	public final double majorPvalue;
	public final double altPvalue;

	public final double majorError;
	public final double altError;

	public final int majorPlusStrandCount;
	public final int majorMinusStrandCount;

	public final int minorPlusStrandCount;
	public final int minorMinusStrandCount;

	public VariantResults(String majorCall, String minorCall, int majorCount,
			int minorCount, int bgCount, double majorPvalue, double altPvalue,
			double majorError, double altError,
			int majorPlusStrandCount, int majorMinusStrandCount, 
			int minorPlusStrandCount, int minorMinusStrandCount) {		
		this.majorCall = majorCall;
		this.minorCall = minorCall;
		
		this.majorCount = majorCount;
		this.minorCount = minorCount;
		this.bgCount = bgCount;
		
		this.majorPvalue = majorPvalue;
		this.altPvalue = altPvalue;

		this.majorError = majorError;
		this.altError = altError;

		this.majorPlusStrandCount = majorPlusStrandCount;
		this.majorMinusStrandCount = majorMinusStrandCount;

		this.minorPlusStrandCount = minorPlusStrandCount;
		this.minorMinusStrandCount = minorMinusStrandCount;
	}

	public double getMajorStrandFreq() {
		if (majorPlusStrandCount == 0 || majorMinusStrandCount == 0) {
			return 0.0;
		}
		if (majorPlusStrandCount < majorMinusStrandCount) {
			return ((double)majorPlusStrandCount) / (majorPlusStrandCount + majorMinusStrandCount);
		}
		return ((double)majorMinusStrandCount) / (majorPlusStrandCount + majorMinusStrandCount);
	}

	public double getMinorStrandFreq() {
		if (minorPlusStrandCount == 0 || minorMinusStrandCount == 0) {
			return 0.0;
		}
		if (minorPlusStrandCount < minorMinusStrandCount) {
			return ((double) minorPlusStrandCount) / (minorPlusStrandCount + minorMinusStrandCount);
		}
		return ((double)minorMinusStrandCount) / (minorPlusStrandCount + minorMinusStrandCount);
	}

}
