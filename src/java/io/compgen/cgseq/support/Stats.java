package io.compgen.cgseq.support;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class Stats {
//	public static final double DEFAULT_DELTA = 1.0e-16;
//	public interface InnerFunc {
//		public double inner(int x);
//	}

	
	/**
	 * http://en.wikipedia.org/wiki/Skellam_distribution
	 * @param k - difference between counts
	 * @param mu1 - expected Poisson mean 1
	 * @param mu2 - expected Poisson mean 2
	 * @return - probability of having the difference (k) between two Poisson counts given distributions means mu1 & mu2
	 */
	
	public static double skellam(int k, double mu1, double mu2) {
		return Math.exp(-mu1-mu2) * Math.pow(mu1/mu2,k/2.0) * BesselI.value(k, 2 * Math.sqrt(mu1 * mu2));
	}

	
	/**
	 * Given two Poisson means, calculate the number where the p-value for the difference (k) is highest.
	 * @param mu1
	 * @param mu2
	 * @return
	 */
	public static int skellamCrossingPoint(double mu1, double mu2) {
		int i = (int) Math.max(mu1, mu2);
		double maxval = 0.0;
		double val = 0.0;

		while (true) {
			val = skellam(i, mu1, mu2);
//			System.err.println("i="+i+", val="+val+ ", maxval="+maxval);
			if (val < maxval) {
				return i+1;
			}
			maxval = val;
			i = i - 1;
		}
	}

	public static double poissonProb(int k, double mu) {
		PoissonDistribution pois = new PoissonDistribution(mu);
		return pois.probability(k);		
	}


	/**
	 * 
	 * @param k - observed successes
	 * @param n - number of trials
	 * @param p - expected rate
	 * @return - one-sided p-value (for two-sided, divide by 2)
	 */
	public static double binomialCumulativeProb(int k, int n, double p) {
		BinomialDistribution binom = new BinomialDistribution(n, p);
		return binom.cumulativeProbability(k);
	}
}
//	/**
//	 * http://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I.CE.B1.2C_K.CE.B1
//	 * http://functions.wolfram.com/Bessel-TypeFunctions/BesselI/02/
//	 * 
//	 * Technically, I think this is the hard way to do this...
//	 * Converting this from FORTRAN would also work: http://www.netlib.org/specfun/ribesl
//	 * (Apache did this for BesselJ)
//	 * 
//	 * @param x
//	 * @param k
//	 * @return
//	 */
//	public static double besselI(final double Z, final int nu) {
//		if (nu < 0) {
//			return besselI(Z, -nu);
//		}
//		
//		/*
//		 * I_v(Z) = sigma(k=0..inf) { (1/(gamma(k + v + 1)* k!)) * (Z/2)^(2k+v) } 
//		 */
//
//		return sumDelta(0, Integer.MAX_VALUE, new InnerFunc() {
//			@Override
//			public double inner(int k) {
//				return (
//						1 / 
//						(Gamma.gamma(k + nu + 1) * CombinatoricsUtils.factorialDouble(k))
//						) 
//						* 
//						Math.pow((Z/2), (2 * k + nu));  
//			}});
//	}
//
//	/**
//	 * Iterate over a range and sum the results (basically an discontinuous integral taken at only integer values)
//	 * This iterates until the return value doesn't change by more than DEFAULT_DELTA
//	 * This is particularly useful when you are iterate to infinity... (or Integer.MAX_VALUE)
//	 * @param from
//	 * @param to
//	 * @param func
//	 * @return
//	 */
//	public static double sumDelta(int from, int to, InnerFunc func) {
//		return sumDelta(from, to, func, DEFAULT_DELTA);
//	}
//	
//	/**
//	 * 
//	 * @param from start value
//	 * @param to end value
//	 * @param func function to use to calculate
//	 * @param delta Stop when the delta between runs is less than this level 
//	 *              Note: the sum needs to change by at least delta once before this 
//	 *              will stop the run. This is required for instances where the "to" is
//	 *              significantly greater than the "from".
//	 * @return
//	 */
//	public static double sumDelta(int from, int to, InnerFunc func, double delta) {
//		double acc = 0.0;
//		double lastacc = acc;
//		
//		boolean hitDelta = false;
//		
//		for (int i=from; i<=to; i++) {
////			System.err.println("i="+i+" acc="+acc+" val=" + func.inner(i));
//			lastacc = acc;
//			acc += func.inner(i);
//
//			if (Double.isNaN(acc)) {
//				return Double.NaN;
//			}
//			
//			double currentDelta = Math.abs(acc - lastacc);
//			
//			if (i > from && currentDelta < delta) {
//				if (hitDelta) {
//					return acc;
//				}
//			} else if (currentDelta > delta) {
////				System.err.println("Hit delta!");
//				hitDelta = true;
//			}
//		}
//		
//		return acc;
//	}
//}
