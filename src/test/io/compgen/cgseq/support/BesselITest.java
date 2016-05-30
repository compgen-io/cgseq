package io.compgen.cgseq.support;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class BesselITest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testValueDoubleDoubleBoolean() {
		innerBesselI(0,1,1.266066);
		innerBesselI(1,1, 0.5651591);
		innerBesselI(2,1, 0.1357477);
		innerBesselI(2,2, 0.6889484);
		innerBesselI(60,45, 405.6376);
		innerBesselI(60,60, 3.320843e12);
		innerBesselI(100,100, 4.641535e21);
		innerBesselI(100,0, 0);
		innerBesselI(100,20, 2.870319e-58);
		innerBesselI(160, 2 * Math.sqrt(160), 1.220186e-108);
		innerBesselI(365, 2 * Math.sqrt(365), 0.0);

		
		innerBesselI(0,1,0.4657596, true);

	}
	
	private void innerBesselI(int order, double x, double expected) {
		innerBesselI(order, x, expected, false, 5);
	}
	private void innerBesselI(int order, double x, double expected, boolean exp) {
		innerBesselI(order, x, expected, exp, 5);
	}
	private void innerBesselI(int order, double x, double expected, boolean exp, int sigfig) {
		double val = BesselI.value(order,x, exp);
//		double old = Stats.besselI(x, order);
		
		System.err.println("BesselI.value("+order+","+x+","+exp+") = " + val + " expected="+expected);// +" old method="+old);
		
		if (Double.toString(val).indexOf('E') > -1) {
			String[] valS = Double.toString(val).split("E");
			String[] expS = Double.toString(expected).split("E");
			
			if (valS.length > 1 || expS.length > 1) {
				assertEquals(valS[1], expS[1]);
			}

			String valS1 = valS[0].substring(0,sigfig+2);
			String expS1 = expS[0].substring(0,sigfig+2);

			
			assertEquals(valS1, expS1);
			
			 
		} else {
			String valS = substr(Double.toString(val), 0,sigfig+2);
			String expS = substr(Double.toString(expected),0,sigfig+2);
			assertEquals(valS, expS);
		}
		
		
	}

	private String substr(String s, int start, int end) {
		if (end > s.length()) {
			end = s.length();
		}
		
		return s.substring(start,  end);
	}
	
	
}
