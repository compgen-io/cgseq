package io.compgen.cgseq.variant;

import io.compgen.cgseq.support.LRUCache;
import io.compgen.cgseq.support.MapCount;
import io.compgen.common.Pair;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCall;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCallOp;

import java.util.List;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class PoissonVariant {
	private final boolean backgroundCorrection;
	private final double expectedAlleleFrequency;
	private final int minQual;
	
	public PoissonVariant(boolean backgroundCorrection, int minQual, double expectedAlleleFrequency) {
		this.backgroundCorrection = backgroundCorrection;
		this.expectedAlleleFrequency = expectedAlleleFrequency;
		this.minQual = minQual;
	}

	public PoissonVariant(boolean backgroundCorrection, int minQual) {
		this.backgroundCorrection = backgroundCorrection;
		this.minQual = minQual;
		this.expectedAlleleFrequency = 0.5;
	}

	public VariantResults calcVariant(List<PileupBaseCall> calls) {
		MapCount<String> counter = new MapCount<String>();
		int rawDepth = 0;
		for (PileupBaseCall call: calls) {
			rawDepth ++;
			if (call.qual > minQual || call.op != PileupBaseCallOp.Match) {
				counter.incr(call.toString());
			}
		}
		
		List<Pair<String, Integer>> sorted = counter.getSortedCounts();

		if (sorted.size()==0) {
			return null;
		}

//		for (Pair<String, Integer> tup: sorted) {
//			System.out.println(tup.one+"\t"+tup.two);
//		}
		
		String majorCall = sorted.get(0).one; 
		String minorCall = ""; 
		
		int major = sorted.get(0).two;
		int minor = 0;
		int bg = 0;
		
		if (sorted.size() > 1) {
			minorCall = sorted.get(1).one;
			minor = sorted.get(1).two;
		}		
		
		if (sorted.size()>2 && backgroundCorrection) {
			bg = sorted.get(2).two;
			major = major - bg;
			minor = minor - bg;
		}

		if (major <= 0 && minor <= 0) {
			return null;
		}

		int majorPlusStrandCount = plusCount(calls, majorCall);
		int majorMinusStrandCount = major - majorPlusStrandCount;
//		int majorStrandCount = majorPlusStrandCount < majorMinusStrandCount ? majorPlusStrandCount: majorMinusStrandCount;

		int minorPlusStrandCount = plusCount(calls, minorCall);
		int minorMinusStrandCount = minor - minorPlusStrandCount;
//		int minorStrandCount = minorPlusStrandCount < minorMinusStrandCount ? minorPlusStrandCount: minorMinusStrandCount;
		
		double majorPvalue = calcPvalue(major, major + minor) + calcPvalue(minor, 1 );
		double majorError = calcCumulativePvalue(minor, (major + minor) * expectedAlleleFrequency) * (1-calcCumulativePvalue(major-1, (major + minor) * expectedAlleleFrequency));
//		double majorStrandPval = calcCumulativePvalue(majorStrandCount, major * 0.5);

		double altPvalue = calcPvalue(minor, (major + minor) * expectedAlleleFrequency) + calcPvalue(major, (major + minor) * expectedAlleleFrequency);
		double altError = calcCumulativePvalue(major, (major + minor)) * calcCumulativePvalue(minor, (major + minor));
//		double minorStrandPval = 1;
		
//		if (minor > 0) {		
//			minorStrandPval = calcCumulativePvalue(minorStrandCount, minor * 0.5);
//		}

		
		// TODO: Also test the quality distributions for major/minor?
		
		return new VariantResults(majorCall, minorCall, rawDepth, majorPvalue - altPvalue);//, majorStrandPval, minorStrandPval);
	}
	
	LRUCache<Pair<Integer, Double>, Double> dCache = new LRUCache<Pair<Integer, Double>, Double>(10000);
	private double calcCumulativePvalue(int observed, double lambda) {
		Pair<Integer, Double> k = new Pair<Integer,Double>(observed, lambda);
		Double pval = dCache.get(k);
		if (pval == null) {
			PoissonDistribution poisson = new PoissonDistribution(lambda);
			pval = poisson.cumulativeProbability(observed);
			dCache.put(k, pval);
		}
		return pval;
	}
	
	LRUCache<Pair<Integer, Double>, Double> pCache = new LRUCache<Pair<Integer, Double>, Double>(10000);
	private double calcPvalue(int observed, double lambda) {
		Pair<Integer, Double> k = new Pair<Integer, Double>(observed, lambda); 
		Double pval = pCache.get(k);
		if (pval == null) {
			PoissonDistribution poisson = new PoissonDistribution(lambda);
			pval = poisson.probability(observed);
			pCache.put(k, pval);
		}
		return pval;
	}
	
	private int plusCount(List<PileupBaseCall> calls, String call) {
		int plus = 0;
		for (PileupBaseCall pbc: calls) {
			if (pbc.qual<minQual) {
				continue;
			}
			if (pbc.toString().equals(call)) {
				if (pbc.plusStrand) {
					plus++;
				}
			}
		}
		return plus;
	}
}
