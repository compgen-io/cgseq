package io.compgen.cgseq.support;

import io.compgen.common.StringUtils;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;

/**
 * This class provides computation methods related to modified Bessel
 * functions of the first kind.
 * 
 * This implementation is based on the ribesl Fortran routine at
 * <a href="http://www.netlib.org/specfun/rjbesl">Netlib</a>. and the BesselJ function
 * from Apache Commons Math (org.apache.commons.math3.special.BesselJ)</p>
 * <p>
 * From the Fortran code: </p>
 * <p>
 * This program is based on a program written by David J. Sookne (2) that
 * computes values of the Bessel functions J or I of real argument and integer
 * order. Modifications include the restriction of the computation to the I
 * Bessel function of non-negative real argument, the extension of the
 * computation to arbitrary positive order, and the elimination of most
 * underflow. An earlier version was published in (3).</p>
 * <p>
 * References:</p>
 * <ul>
 * <li>"A Note on Backward Recurrence Algorithms," Olver, F. W. J., and Sookne,
 * D. J., Math. Comp. 26, 1972, pp 941-947.</li>
 * <li>"Bessel Functions of Real Argument and Integer Order," Sookne, D. J., NBS
 * Jour. of Res. B. 77B, 1973, pp 125-132.</li>
 * <li>"ALGORITHM 597, Sequence of Modified Bessel Functions of the First Kind," 
 * Cody, W. J., Trans. Math. Soft.,1983, pp. 242-245.</li>
 * </ul> </p>
 * @since 3.4
 */

public class BesselI implements UnivariateFunction {

	// Machine dependent constants 

    /**
     * 10.0^K, where K is the largest integer such that ENTEN is
     * machine-representable in working precision
     */
    private static final double ENTEN = 1.0e308;

    /**
     * Decimal significance desired. Should be set to (INT(log_{10}(2) * (it)+1)).
     * Setting NSIG lower will result in decreased accuracy while setting
     * NSIG higher will increase CPU time without increasing accuracy.
     * The truncation error is limited to a relative error of
     * T=.5(10^(-NSIG)).
     */
    private static final double ENSIG = 1.0e16;

    /** 10.0 ** (-K) for the smallest integer K such that K >= NSIG/4 */
    private static final double RTNSIG = 1.0e-4;

    /** Smallest ABS(X) such that X/4 does not underflow */
    private static final double ENMTEN = 8.90e-308;

    /**
     * Upper limit on the magnitude of x. If abs(x) = n, then at least
     * n iterations of the backward recursion will be executed. The value of
     * 10.0 ** 4 is used on every machine.
     */
    private static final double XLARGE = 1.0e4;

    /**
     * Largest working precision argument that the library
     * EXP routine can handle and upper limit on the
     * magnitude of X when IZE=1; approximately 
     * LOG(beta**maxexp)
     */
    private static final double EXPARG = 709;

    /**
	 * Decimal significance desired.  Should be set to
	 * INT(LOG10(2)*it+1).  Setting NSIG lower will result
	 * in decreased accuracy while setting NSIG higher will
	 * increase CPU time without increasing accuracy.  The
	 * truncation error is limited to a relative error of
	 * T=.5*10**(-NSIG).     
	 */
    
    private static final double NSIG = 16;

    private static final double CONST = 1.585;

    /** Order of the function computed when {@link #value(double)} is used */
    private final double order;

    /**
     * Create a new BesselI with the given order.
     *
     * @param order order of the function computed when using {@link #value(double)}.
     * @param expScale exponentially scaled I's are to be calculated (var "ize").
     */
    public BesselI(double order) {
        this.order = order;
    }

	@Override
	public double value(double x) {
		return BesselI.value(order, x);
	}

	public static double value(double order, double x) {
		return value(order, x, false);
	}

	/**
	 * Returns modified first Bessel function \(I_{order}(x)\)
	 * @param order
	 * @param x
     * @return Value of the Bessel function of the modified first kind, \(I_{order}(x)\)
     * @throws MathIllegalArgumentException if {@code x} is too large relative to {@code order}
     * @throws ConvergenceException if the algorithm fails to converge
	 */
	public static double value(double order, double x, boolean expScale) {
		final int n = (int) order;
		final double alpha = order - n;
		final int nb = n + 1;
		
		final BesselIResult res = riBesl(x, alpha, nb, expScale);
		if (res == null) {
			return -1;
		}
		
		// res.vals is 1-based, so we'll keep that here.
//		System.err.println(res.nVals+"/"+StringUtils.join(";", res.vals));
        if (res.nVals >= nb) {
            return res.vals[n+1];
        } else if (res.nVals < 0) {
            throw new MathIllegalArgumentException(LocalizedFormats.BESSEL_FUNCTION_BAD_ARGUMENT,order, x);
        } else if (FastMath.abs(res.vals[res.nVals]) < 1e-100) {
            return res.vals[n+1]; // underflow; return value (will be zero)
        }
        throw new ConvergenceException(LocalizedFormats.BESSEL_FUNCTION_FAILED_CONVERGENCE, order, x);
	}
	
	
	public static class BesselIResult {
		private final double[] vals;
		private final int nVals;

		public BesselIResult(double[] b, int n) {
			vals = MathArrays.copyOf(b, b.length);
			nVals = n;
		}
		
		public double[] getVals() {
			return MathArrays.copyOf(vals, vals.length);
		}
		
		public int getnVals() {
			return nVals;
		}
		
	}	
	
	public static BesselIResult riBesl(final double x, final double alpha, final int nb, final boolean ize) {
		final double[] b = new double[nb+1]; // fortran is 1-based, so we'll keep that convention here.
		
		int ncalc = 0;
		
		// fail fast
		if (nb <= 0 || x < 0 || alpha < 0 || alpha >= 1 || (!ize && x > EXPARG) || (ize && x > XLARGE)) {
			ncalc = Math.min(nb, 0)-1;
			return new BesselIResult(MathArrays.copyOf(b, b.length), ncalc);
		}

		ncalc = nb;
		final int magx = (int) x;

		// initialize
        for (int i = 0; i <= nb; ++i) {
            b[i] = 0;
        }

		// Use 2-term ascending series for small X ?? This is the large X branch...

        if (x >= RTNSIG) {
//        	System.err.println("x >= RTNSIG =>" + x + " >= " + RTNSIG);

			// Initialize the forward sweep, the P-sequence of Olver
			int nbmx = nb-magx;
			int n = magx + 1;
			double en = (n + n) + (alpha + alpha);
			double plast = 1;
			double p = en / x;
			
			// Calculate general significance test
			double test = ENSIG + ENSIG;
			if ((2 * magx) > (5 * NSIG)) {
				test = FastMath.sqrt(test * p);
			} else {
				test = test / FastMath.pow(CONST, magx);
			}
			
			// line 228
			
			
			
			boolean goto120 = false;
			if (nbmx >= 3) {
//				System.err.println("nbmx >= 3 ; "+ nbmx +" >= " + 3);
				// Calculate P-sequence until N = NB-1; check for possible overflow
				double TOVER = ENTEN / ENSIG;
				int nstart = magx + 2;
				int nend = nb - 1;
				
				for (int k = nstart; k<=nend; k++) { // to 100
					n = k;
					en = en + 2;
					double pold = plast;
					plast = p;
					p = en * plast / x + pold;
					
					if (p > TOVER) {
						// to avoid overflow, divide p-sequence by tover.
						// p-sequence until abs(p) > 1;
						
						TOVER = ENTEN;
						p = p / TOVER;
						plast = plast / TOVER;
						double psave = p;
						double psavel = plast;
						nstart = n + 1; 

						// 60
						do {
							n = n + 1;
							en = en + 2;
							pold = plast;
							plast = p;
							p = en * plast / x + pold;
						} while (p <= 1.0);
						
						double tempb = en / x;
						
						// Calculate backward test, and find NCALC, the highest N
						// such that the test is passed
						
						test = pold * plast / ENSIG;
						test = test * (0.5 - 0.5 / (tempb * tempb));
						p = plast * TOVER;
						n = n - 1;
						en = en - 2;
						nend = FastMath.min(nb, n);
						
						for (int l = nstart; l <= nend && psave * psavel <= test; l++) {
							ncalc = l;
							pold = psavel;
							psavel = psave;
							psave = en * psavel / x + pold;
						}

						if (psave * psavel > test) {
							// 90
							ncalc = ncalc - 1;
						} else {
							// 80
							ncalc = nend + 1;
						}
						goto120 = true;
					} // end if p > tover 
				} // 100 - end do loop
					
				if (!goto120) {
					n = nend;
					en = (double) n+n + (alpha + alpha);
					// Calculate special significance test for nbmx > 2
					test = FastMath.max(test, FastMath.sqrt(plast * ENSIG) * FastMath.sqrt(p+p));
				}
			} else {
//				System.err.println("nbmx < 3 ; "+ nbmx +" < " + 3);
				// endif nbmx >= 3
			}
			
			// 110 
			// Calculate P-sequence until significance test passed
			
			if (!goto120) {
				double pold;
				do {
//					System.err.println("p = "+p+ "; test = "+ test);
					n = n + 1;
					en = en + 2;
					pold = plast;
					plast = p;
					p = en * plast / x + pold;
				} while (p < test);
			}
				
			// 120
			
			n = n + 1;
			en = en + 2.0;
			double tempb = 0.0;
			double tempa = 1.0 / p;
			double em = n - 1.0;
			double empal = em + alpha;
			double emp2al = (em - 1) + (alpha + alpha);
			double sum = tempa * empal * emp2al / em;
			int nend = n - nb;
			
			boolean goto230 = false;
			
			// line 309
			if (nend < 0) {
//				System.err.println("nend < 0; nend="+ nend);
				
				// N .LT. NB, so store B(N) and set higher orders to zero.

				b[n] = tempa;
				nend = -nend;
				for (int i = 1; i <= nend; i++) {
					b[n+i] = 0;
				}
			} else {
				if (nend > 0) {
//					System.err.println("nend > 0; nend="+ nend);
					// Recur backward via difference equation, calculating (but
					// not storing) B(N), until N = NB.
						
					for (int i = 1; i <=nend; i++) {
						n = n - 1;
						en = en - 2;
						double tempc = tempb;
						tempb = tempa;
						tempa = (en * tempb) / x + tempc;
						em = em - 1;
						emp2al = emp2al - 1;
						if (n == 1) {
							break; // goto 150
						}
						if (n == 2) {
							emp2al = 1;
						}
						empal = empal - 1;
						sum = (sum + tempa*empal) * emp2al / em;
					}
				}
				
				// Store B(NB)
				// 150
				b[n] = tempa;
				if (nb <= 1) {
					sum = (sum + sum) + tempa;
					// goto 230;
					goto230 = true;
				} else {
					// Calculate and Store B(NB-1)
					n = n - 1;
					en = en - 2;
					b[n] = en * tempa / x + tempb;
					if (n == 1) {
						// goto 220
						sum = (sum + sum) + b[1]; // copy of line 220, then skip to 230
						goto230 = true;
					} else {
						em = em - 1;
						emp2al = emp2al - 1;
						if (n == 2) {
							emp2al = 1;
						}
						
						empal = empal - 1;
						sum = (sum + b[n] * empal) * emp2al / em;
					}
				}
				
				//
				
				if (!goto230) {
					nend = n - 2;
					if (nend > 0) {
						// Calculate via difference equation and store B(N), until N = 2.
						
						for (int i=1; i<=nend; i++) {
							n = n - 1;
							en = en - 2;
							b[n] = (en * b[n+1]) / x + b[n+2];
							em = em -1;
							emp2al = emp2al - 1;
							if (n == 2) {
								emp2al = 1;
							}
							empal = empal - 1;
							sum = (sum + b[n] * empal) * emp2al / em;
						} // 200
					}
					b[1] = 2 * empal * b[2] / x + b[3];
					sum = (sum + sum) + b[1];
				}
				
				// 230
				if (alpha != 0.0) {
					sum = sum * Gamma.gamma(1+alpha) * FastMath.pow((x*0.5), -alpha);
				}
				
				if (!ize) { // ize == 1
					sum = sum * FastMath.exp(-x);
				}
				
				tempa = ENMTEN;
				if (sum > 1) {
					tempa = tempa * sum;
				}
				
				for (int i=1; i<=nb; i++) {
					if (b[i] < tempa) {
						b[i] = 0;
					}
					b[i] = b[i] / sum;
				} // 260
				// done with this branch - return b[]
			}
		} else {
//        	System.err.println("x < RTNSIG =>" + x + " < " + RTNSIG);
			// Two-term ascending series for small X.
			// Line ~ 262 (actually 395)
			double tempa = 1.0;
			double empal = 1.0 + alpha;
			double halfx = 0.0;
			
			if (x > ENMTEN) {
				halfx = 0.5 * x;
			}
			
			if (alpha != 0.0) {
				tempa = FastMath.pow(halfx, alpha) / Gamma.gamma(empal);
			}
			
			if (ize) {
				tempa = tempa * FastMath.exp(-x);
			}
			
			double tempb = 0.0;
			
			if ((x + 1.0) > 1.0) {
				tempb = halfx * halfx;
			}
			
			b[1] = tempa + tempa * tempb / empal;
			
			if (x != 0.0 && b[0] == 0.0) {
				ncalc = 0;
			}
			
			if (nb > 1) {
				if (x > 0.0) {
					// calculate higher-order functions (loop-310 removed since we already initialize b[] to zeros)
					double tempc = halfx;
					double tover = (ENMTEN + ENMTEN) / x;
					if (tempb != 0.0) {
						tover = ENMTEN / tempb;
					}
					
					for (int i=2; i<=nb; i++) {
						tempa = tempa / empal;
						empal = empal + 1.0;
						tempa = tempa * tempc;
						if (tempa <= (tover * empal)) { 
							tempa = 0.0;
						}
						b[i] = tempa + tempa*tempb / empal;
						if (b[i] == 0.0 && ncalc > i) {
							ncalc = i-1;
						}
					}
				}
			}
		} 
		return new BesselIResult(MathArrays.copyOf(b, b.length), ncalc);

	}

}
