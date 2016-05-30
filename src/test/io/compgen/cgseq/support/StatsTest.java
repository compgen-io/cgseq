package io.compgen.cgseq.support;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class StatsTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testSkellamCrossingPoint() {
		System.err.println("Skellam CP: (15, 30) => k=" + Stats.skellamCrossingPoint(15,  30) + ", val="+Stats.skellam(Stats.skellamCrossingPoint(15,  30), 15, 30));
		System.err.println("Skellam CP: (30, 15) => k=" + Stats.skellamCrossingPoint(30,  15) + ", val="+Stats.skellam(Stats.skellamCrossingPoint(30, 15), 30, 15));
		System.err.println("Skellam CP: (60, 30) => k=" + Stats.skellamCrossingPoint(60,  30) + ", val="+Stats.skellam(Stats.skellamCrossingPoint(60, 30), 60, 30));
	}
	
	@Test
	public void testSkellam() {
		System.err.println("k=4, mu1=1.0, mu2=2.0 => "+ Stats.skellam(4, 1.0, 2.0));
		System.err.println("k=3, mu1=1.0, mu2=2.0 => "+ Stats.skellam(3, 1.0, 2.0));
		System.err.println("k=2, mu1=1.0, mu2=2.0 => "+ Stats.skellam(2, 1.0, 2.0));
		System.err.println("k=1, mu1=1.0, mu2=2.0 => "+ Stats.skellam(1, 1.0, 2.0));
		System.err.println("k=0, mu1=1.0, mu2=2.0 => "+ Stats.skellam(0, 1.0, 2.0));
//		System.err.println("k=-1, mu1=1.0, mu2=2.0 => "+ Stats.skellam(-1, 1.0, 2.0));
//		System.err.println("k=-2, mu1=1.0, mu2=2.0 => "+ Stats.skellam(-2, 1.0, 2.0));
//		System.err.println("k=-3, mu1=1.0, mu2=2.0 => "+ Stats.skellam(-3, 1.0, 2.0));
//		System.err.println("k=-4, mu1=1.0, mu2=2.0 => "+ Stats.skellam(-4, 1.0, 2.0));
		// These vals come from R's calculation...
		assertEquals(0.003056575, Stats.skellam(4, 1.0, 2.0), 0.0000001);
		assertEquals(0.01337568, Stats.skellam(3, 1.0, 2.0), 0.0000001);
		assertEquals(0.04624018, Stats.skellam(2, 1.0, 2.0), 0.0000001);
		assertEquals(0.1192317, Stats.skellam(1, 1.0, 2.0), 0.0000001);
		assertEquals(0.2117121, Stats.skellam(0, 1.0, 2.0), 0.0000001);
//		assertEquals(0.2384634, Stats.skellam(-1, 1.0, 2.0), 0.0000001);
//		assertEquals(0.1849607, Stats.skellam(-2, 1.0, 2.0), 0.0000001);
//		assertEquals(0.1070054, Stats.skellam(-3, 1.0, 2.0), 0.0000001);
//		assertEquals(0.0489052, Stats.skellam(-4, 1.0, 2.0), 0.0000001);

		System.err.println("k=160, mu1=160, mu2=1 => " + Stats.skellam(160, 160.0, 1));
		
	}
	
	@Test
	public void testBinomial() {
		System.err.println("25,100,0.5 => " + Stats.binomialCumulativeProb(25, 100, 0.5));
	}

//	@Test
//	public void testBesselI() {
//		assertEquals(0.5651591, Stats.besselI(1,1), 0.000001);
//		assertEquals(1.590637, Stats.besselI(2,1), 0.000001);
//		assertEquals(9.759465, Stats.besselI(4,1), 0.000001);
//	}
//	
//	@Test
//	public void testSumIntIntInnerFunc() {
//		assertEquals(2 + 4 + 6 + 8 + 10 + 12 + 14 + 16 + 18 + 20, Stats.sumDelta(1, 10, new InnerFunc(){
//			@Override
//			public double inner(int x) {
//				return x * 2;
//			}}), 0);
//		
//	}
}
