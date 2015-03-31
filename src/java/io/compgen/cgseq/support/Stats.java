package io.compgen.cgseq.support;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;

public class Stats {
	public static final double DEFAULT_DELTA = 1.0e-16;
	public interface InnerFunc {
		public double inner(int x);
	}

	/**
	 * http://en.wikipedia.org/wiki/Skellam_distribution
	 * @param k
	 * @param mu1
	 * @param mu2
	 * @return
	 */
	
	public static double skellam(int k, double mu1, double mu2) {
		return Math.exp(-mu1-mu2) * Math.pow(mu1/mu2,k/2.0) * besselI(2 * Math.sqrt(mu1 * mu2), k);
	}

	public static int skellamCrossingPoint(double mu1, double mu2) {
		int i = (int) Math.max(mu1, mu2);
		double maxval = 0.0;
		double val = 0.0;

		while (true) {
			val = skellam(i, mu1, mu2);
			System.err.println("i="+i+", val="+val+ ", maxval="+maxval);
			if (val < maxval) {
				return i+1;
			}
			maxval = val;
			i = i - 1;
		}
	}
	
	/**
	 * http://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions_:_I.CE.B1.2C_K.CE.B1
	 * http://functions.wolfram.com/Bessel-TypeFunctions/BesselI/02/
	 * 
	 * Technically, I think this is the hard way to do this...
	 * Converting this from FORTRAN would also work: http://www.netlib.org/specfun/ribesl
	 * (Apache did this for BesselJ)
	 * 
	 * @param x
	 * @param k
	 * @return
	 */
	public static double besselI(final double Z, final int nu) {
		if (nu < 0) {
			return besselI(Z, -nu);
		}
		
		/*
		 * I_v(Z) = sigma(k=0..inf) { (1/(gamma(k + v + 1)* k!)) * (Z/2)^(2k+v) } 
		 */

		return sum(0, Integer.MAX_VALUE, new InnerFunc() {
			@Override
			public double inner(int k) {
				return (
						1 / 
						(Gamma.gamma(k + nu + 1) * CombinatoricsUtils.factorialDouble(k))
						) 
						* 
						Math.pow((Z/2), (2 * k + nu));  
			}});
	}

	/**
	 * Iterate over a range and sum the results
	 * @param from
	 * @param to
	 * @param func
	 * @return
	 */
	public static double sum(int from, int to, InnerFunc func) {
		return sum(from, to, func, DEFAULT_DELTA);
	}
	
	public static double sum(int from, int to, InnerFunc func, double delta) {
		double acc = 0.0;
		double lastacc = acc;
		for (int i=from; i<=to; i++) {
			lastacc = acc;
			acc += func.inner(i);

			if (acc == Double.NaN) {
				return Double.NaN;
			}
			
			if (i > from && Math.abs(acc - lastacc) < delta) {
				return acc;
			}
		}
		return acc;
	}
}
