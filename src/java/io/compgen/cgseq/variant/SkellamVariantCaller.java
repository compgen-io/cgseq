package io.compgen.cgseq.variant;

import io.compgen.cgseq.support.LRUCache;
import io.compgen.cgseq.support.MapCount;
import io.compgen.cgseq.support.Stats;
import io.compgen.common.ListBuilder;
import io.compgen.common.Pair;
import io.compgen.common.StringUtils;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCall;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCallOp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

public class SkellamVariantCaller implements VariantCaller {
	
	private class SkellamMemoKey {
		public final int k;
		public final double mu1;
		public final double mu2;
		
		private SkellamMemoKey(int k, double mu1, double mu2)
		{
			this.k = k;
			this.mu1 = mu1;
			this.mu2 = mu2;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + k;
			long temp;
			temp = Double.doubleToLongBits(mu1);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(mu2);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			SkellamMemoKey other = (SkellamMemoKey) obj;
			if (!getOuterType().equals(other.getOuterType())) {
				return false;
			}
			if (k != other.k) {
				return false;
			}
			if (Double.doubleToLongBits(mu1) != Double
					.doubleToLongBits(other.mu1)) {
				return false;
			}
			if (Double.doubleToLongBits(mu2) != Double
					.doubleToLongBits(other.mu2)) {
				return false;
			}
			return true;
		}

		private VariantCaller getOuterType() {
			return SkellamVariantCaller.this;
		}
	}
	
	private final MannWhitneyUTest mwut = new MannWhitneyUTest();

	private final boolean backgroundCorrection;
	private final double expectedHeterozygousFrequency;
	private final double expectedHomozygousFrequency;
	private final int minQual;
	private final int minDepth;
	
	private LRUCache<SkellamMemoKey, Double> cache = new LRUCache<SkellamMemoKey, Double>(10000);
	private List<String> infoFields = Collections.unmodifiableList(new ListBuilder<String>()
										.add("INDEL")
										.add("DPR")
										.list());
	
	private List<String> formatFields = Collections.unmodifiableList(new ListBuilder<String>()
										.add("GT")
										.add("DP")
										.add("DP4")
										.add("DV")
										.add("BG")
										.add("SB")
										.add("RPB")
										.add("MSF")
										.add("DEBUG")
										.list());
	
	public SkellamVariantCaller(boolean backgroundCorrection, int minQual, int minDepth, double expectedAlleleFrequency, double expectedHomozygousFrequency) {
		this.backgroundCorrection = backgroundCorrection;
		this.expectedHomozygousFrequency = expectedHomozygousFrequency;
		this.expectedHeterozygousFrequency = expectedAlleleFrequency;
		this.minQual = minQual;
		this.minDepth = minDepth;
	}

	public SkellamVariantCaller(boolean backgroundCorrection, int minQual, int minDepth) {
		this(backgroundCorrection, minQual, minDepth, 0.5, 1.0);
	}

	/* (non-Javadoc)
	 * @see io.compgen.cgseq.variant.VariantCaller#calcVariant(java.util.List)
	 */
	@Override
	public VariantResults calcVariant(List<PileupBaseCall> calls, String ref) {
		MapCount<String> counter = new MapCount<String>();
		int rawDepth = 0;
		
		for (PileupBaseCall call: calls) {
			rawDepth++;
			
			// mpileup doesn't report out indel quality scores, so we just accept them all.
			if (call.qual > minQual || call.op == PileupBaseCallOp.Ins || call.op == PileupBaseCallOp.Del) {
				counter.incr(call.toString());
			}
		}
		
		List<Pair<String, Integer>> sorted = counter.getSortedCounts();

		if (rawDepth < minDepth) {
			return null;
		}
		
		String majorCall = sorted.get(0).one; 
		String minorCall = "";
		
		int major = sorted.get(0).two;
		int minor = 0;
		int bg = 0;
		
		if (sorted.size() > 1) {
			minorCall = sorted.get(1).one;
			minor = sorted.get(1).two;
		}		
		
		if (sorted.size() > 2 && backgroundCorrection) {
			bg = sorted.get(2).two;
			major = major - bg;
			minor = minor - bg;
		}

		if (major <= 0 && minor <= 0) {
			return null;
		}

		int majorPlusStrandCount = plusCount(calls, majorCall);
		int majorMinusStrandCount = major - majorPlusStrandCount;
		double majorMSF = (double) Math.min(majorPlusStrandCount, majorMinusStrandCount) / major;
		double majorSB = getBinomialProb(Math.min(majorPlusStrandCount, majorMinusStrandCount), major, 0.5);
		
		int minorPlusStrandCount = 0;
		int minorMinusStrandCount = 0;
		double minorMSF = 0.0;
		double minorSB = 0.0;
		
		if (minor > 0) {
			minorPlusStrandCount = plusCount(calls, minorCall);
			minorMinusStrandCount = minor - minorPlusStrandCount;
			minorMSF = (double) Math.min(minorPlusStrandCount, minorMinusStrandCount) / minor;
			minorSB = getBinomialProb(Math.min(minorPlusStrandCount, minorMinusStrandCount), minor, 0.5);
		}

		int diff = major - minor;
		double hom = (major + minor) * expectedHomozygousFrequency;
		double het = (major + minor) * expectedHeterozygousFrequency;

		double probHom;
		double probHet;
		
		double rpb = 0.0;
		
		boolean pois = false; // did we use a Poisson test or the Skellam test.
		
		try {
			probHom = getSkellamProb(diff, hom, 1); // probability of hom call assuming 1 alt-call (seq error).
			probHet = getSkellamProb(diff, het, het);
		} catch (MathIllegalArgumentException ex) {
			// if the counts are too high, revert to a plain Poisson based test
			
			// probHom = Poisson(major; expected hom count)
			// probHet = Poisson(minor; expected het count)

			pois = true;
			probHom = getPoissonProb(major, hom);
			probHet = getPoissonProb(minor, het);
		}
		
		if (probHom == 0.0 || probHet == 0.0) {
			// if probHom or probHet is 0.0, that is an underflow in the BesselI function
			// This means the difference between hom and het calls is great and the depth
			// is also high.
			// This should only happen for high depth bases (~300X) YMMV. And at this size, the 
			// Skellam dist isn't as critical

			// reverting to a plain Poisson based test

			pois = true;
			probHom = getPoissonProb(major, hom);
			probHet = getPoissonProb(minor, het);
			
		}
		
		
		VariantResults results;
		
		boolean isHet = false;
		
		if (probHet >= probHom && minor > 0) {
			results = new VariantResults(majorCall, minorCall, rawDepth, probHom);
			isHet = true;
			if (majorCall.length() > 1 || minorCall.length()>1) {
				results.addInfo("INDEL");
			}

			List<Integer> majorPos = getReadPos(calls, majorCall);
			double[] majpos = new double[majorPos.size()];
			for (int i=0; i<majpos.length; i++){
				majpos[i] = majorPos.get(i);
			}

			List<Integer> minorPos = getReadPos(calls, minorCall);
			double[] minpos = new double[minorPos.size()];
			for (int i=0; i<minpos.length; i++){
				minpos[i] = minorPos.get(i);
			}
			
			rpb = mwut.mannWhitneyUTest(majpos, minpos);
//			System.err.println("MWUT => "+StringUtils.join(",", majpos)+ " vs " + StringUtils.join(",", minpos) +" U= " + mwut.mannWhitneyU(majpos, minpos) +", p=" + rpb);
			
		} else {
			results = new VariantResults(majorCall, null, rawDepth, probHet);
			if (majorCall.length() > 1) {
				results.addInfo("INDEL");
			}
		} 

		results.addFormat("BG", bg);
		
		// high quality depth
		results.addFormat("DP", rawDepth);
		
		results.addFormat("DEBUG", majorCall+","+minorCall+","+ref+","+probHom+","+probHet+","+major+","+minor+(pois? ",pois":",skel"));
		
		// high quality non-reference bases
		if (majorCall.equals(ref)) {
			// major call is ref, minor call is alt
 			results.addFormat("DV", minor);
 			if (isHet) {
 				results.addInfo("DPR", major+","+minor);
 				results.addFormat("GT", "0/1");
 				results.addFormat("DP4", ""+majorPlusStrandCount+","+majorMinusStrandCount+","+minorPlusStrandCount+","+minorMinusStrandCount);
 				results.addFormat("MSF", majorMSF + "," + minorMSF);
 				results.addFormat("SB", majorSB +"," + minorSB);
 				results.addFormat("RPB", rpb);
 			} else {
				results.addInfo("DPR", major);
 				results.addFormat("GT", "0/0");
				results.addFormat("DP4", ""+majorPlusStrandCount+","+majorMinusStrandCount);
 				results.addFormat("MSF", majorMSF);
 				results.addFormat("SB", majorSB);
			}
		} else if (minorCall.equals(ref) || minorCall.equals("")) {
			// minor call is ref, major call is alt
			results.addFormat("DV", major);

			// even though this isn't a het, we have to return the "minor" values for the ref allele
			results.addInfo("DPR", minor+","+major);
			results.addFormat("MSF", minorMSF + "," + majorMSF);
			results.addFormat("DP4", ""+minorPlusStrandCount+","+minorMinusStrandCount+","+majorPlusStrandCount+","+majorMinusStrandCount);
			results.addFormat("SB", minorSB +"," + majorSB);

			if (isHet) {
 				results.addFormat("GT", "0/1");
 				results.addFormat("RPB", rpb);
 			} else {
 				results.addFormat("GT", "1/1");
			}
		} else {
			// het with two alt alleles
			results.addInfo("DPR", "0,"+major+","+minor);
			results.addFormat("GT", "1/2");
			results.addFormat("DV", major+minor);
			results.addFormat("DP4", "0,0,"+majorPlusStrandCount+","+majorMinusStrandCount+","+minorPlusStrandCount+","+minorMinusStrandCount);
			results.addFormat("MSF", "0,"+majorMSF+","+minorMSF);
			results.addFormat("SB", "0,"+majorSB +"," + minorSB);
			results.addFormat("RPB", rpb);
		}
		
		
		
//		##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
//		##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">

//		##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
//		##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
//		##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
//		##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
//		##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
//		##INFO=<ID=DPR,Number=R,Type=Integer,Description="Number of high-quality bases observed for each allele">

//		##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
//		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
//		##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of high-quality non-reference bases">
//		##FORMAT=<ID=DPR,Number=R,Type=Integer,Description="Number of high-quality bases observed for each allele">
//		##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases">
//		##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">

		
		return results;

	}

	/**
	 * 
	 * @param k - count
	 * @param mu - expected mean
	 * @return
	 */
	private double getPoissonProb(int k, double mu) {
		return Stats.poissonProb(k, mu);
	}

	/**
	 * 
	 * @param x - observed successes
	 * @param n - number of trials
	 * @param p - expected rate
	 * @return
	 */
	private double getBinomialProb(int k, int n, double p) {
		return Stats.binomialCumulativeProb(k, n, p);
	}

	
	/**
	 * 
	 * @param k - observed difference between two Poisson values
	 * @param mu1 - mean 1
	 * @param mu2 - mean 2
	 * @return
	 */
	private double getSkellamProb(int k, double mu1, double mu2) {
		SkellamMemoKey memo = new SkellamMemoKey(k, mu1, mu2);
		if (cache.containsKey(memo)) {
			return cache.get(memo);
		}

		double val = Stats.skellam(k, mu1, mu2);
		if (Double.isNaN(val)) {
			val = 0.0;
		}
		cache.put(memo, val);
		
		return val;
	}

	private int plusCount(List<PileupBaseCall> calls, String call) {
		int plus = 0;
		for (PileupBaseCall pbc: calls) {
			if (pbc.qual<minQual && pbc.op == PileupBaseCallOp.Match) {
				continue;
			}
			if (pbc.matches(call)) {
				if (pbc.plusStrand) {
					plus++;
				}
			}
		}
		return plus;
	}

	private List<Integer> getReadPos(List<PileupBaseCall> calls, String call) {
		List<Integer> readpos = new ArrayList<Integer>();
		for (PileupBaseCall pbc: calls) {
			if (pbc.qual<minQual && pbc.op == PileupBaseCallOp.Match) {
				continue;
			}
			if (pbc.matches(call)) {
				readpos.add(pbc.readPos);
			}
		}
		return readpos;
	}

	@Override
	public List<String> getInfoFields() {
		return infoFields;
	}

	@Override
	public List<String> getFormatFields() {
		return formatFields;
	}

	@Override
	public String getInfoFieldDescription(String k) {
		switch(k) {
		case "INDEL":
			return "Number=0,Type=Flag,Description=\"Variant is an IN-DEL.\">";
		case "DPR":
			return "Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\"";
		}
		return null;
	}

	@Override
	public String getFormatFieldDescription(String k) {
		switch(k) {
		case "GT":
			return "Number=1,Type=String,Description=\"Genotype calls\">";
		case "DP4":
			return "Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases (or alt1-fwd, alt1-rev, alt2-fwd, alt2-rev)\">";
		case "DP":
			return "Number=1,Type=Integer,Description=\"# high-quality bases (raw-depth)\"";
		case "BG":
			return "Number=1,Type=Integer,Description=\"Background calls\"";
		case "DV":
			return "Number=1,Type=Integer,Description=\"# high-quality non-reference bases\"";
		case "SB":
			return "Number=R,Type=Integer,Description=\"Strand-bias p-value for each allele (Binomial)\"";
		case "RPB":
			return "Number=R,Type=Integer,Description=\"Read position bias p-value (Mann-Whitney U Test)\"";
		case "MSF":
			return "Number=R,Type=Integer,Description=\"Minor strand frequency\"";
		}
		return null;
	}
}
