package io.compgen.cgseq.variant;

import io.compgen.cgseq.support.LRUCache;
import io.compgen.cgseq.support.MapCount;
import io.compgen.cgseq.support.Stats;
import io.compgen.common.MapBuilder;
import io.compgen.common.Pair;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCall;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCallOp;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

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
	
	private final boolean backgroundCorrection;
	private final double expectedHeterozygousFrequency;
	private final double expectedHomozygousFrequency;
	private final int minQual;
	private final int minDepth;
	
	private LRUCache<SkellamMemoKey, Double> cache = new LRUCache<SkellamMemoKey, Double>(10000);
	
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
			if (call.qual > minQual || call.op != PileupBaseCallOp.Match) {
				counter.incr(call.toString());
			}
		}
		
		List<Pair<String, Integer>> sorted = counter.getSortedCounts();

		if (sorted.size() < minDepth) {
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
//		int majorStrandCount = majorPlusStrandCount < majorMinusStrandCount ? majorPlusStrandCount: majorMinusStrandCount;

		int minorPlusStrandCount = plusCount(calls, minorCall);
		int minorMinusStrandCount = minor - minorPlusStrandCount;
//		int minorStrandCount = minorPlusStrandCount < minorMinusStrandCount ? minorPlusStrandCount: minorMinusStrandCount;

		
		int diff = major - minor;
		double hom = (major + minor) * expectedHomozygousFrequency;
		double het = (major + minor) * expectedHeterozygousFrequency;
		
		double probHom = getSkellamProb(diff, hom, 1);
		double probHet = getSkellamProb(diff, het, het);
		
		VariantResults results;
		if (probHom > probHet) {
			results = new VariantResults(majorCall, null, rawDepth, probHom);
		} else {
			results = new VariantResults(majorCall, minorCall, rawDepth, probHet);
		}
		
		// high quality depth
		results.addFormat("DP", ""+major+minor);
		
		// high quality non-reference bases
		if (majorCall.equals(ref)) {
			results.addFormat("DV", ""+minor);
			results.addFormat("DP4", ""+majorPlusStrandCount+","+majorMinusStrandCount+","+minorPlusStrandCount+","+minorMinusStrandCount);
			results.addInfo("DPR", major+","+minor);
		} else if (minorCall.equals(ref)) {
			results.addFormat("DV", ""+major);
			results.addFormat("DP4", ""+minorPlusStrandCount+","+minorMinusStrandCount+","+majorPlusStrandCount+","+majorMinusStrandCount);
			results.addInfo("DPR", minor+","+major);
		} else {
			results.addFormat("DV", ""+major+minor);
			results.addFormat("DP4", "0,0,"+(majorPlusStrandCount+minorPlusStrandCount)+","+(majorMinusStrandCount+minorMinusStrandCount));
			results.addInfo("DPR", "0,"+major+","+minor);
		}
		
		
//		##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
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
	
	private double getSkellamProb(int k, double mu1, double mu2) {
		SkellamMemoKey memo = new SkellamMemoKey(k, mu1, mu2);
		if (cache.containsKey(memo)) {
			return cache.get(memo);
		}

		double val = Stats.skellam(k, mu1, mu2);
		cache.put(memo, val);
		
		return val;
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

	private List<Integer> getReadPosStrand(List<PileupBaseCall> calls, String call, boolean plusStrand) {
		List<Integer> readpos = new ArrayList<Integer>();
		for (PileupBaseCall pbc: calls) {
			if (pbc.qual<minQual) {
				continue;
			}
			if (pbc.toString().equals(call)) {
				if (pbc.plusStrand == plusStrand) {
					readpos.add(pbc.readPos);
				}
			}
		}
		return readpos;
	}

	private List<Integer> getReadPos(List<PileupBaseCall> calls, String call) {
		List<Integer> readpos = new ArrayList<Integer>();
		for (PileupBaseCall pbc: calls) {
			if (pbc.qual<minQual) {
				continue;
			}
			if (pbc.toString().equals(call)) {
				readpos.add(pbc.readPos);
			}
		}
		return readpos;
	}

	@Override
	public Map<String, String> getInfoFields() {
		return new MapBuilder<String, String>()
				.put("DPR", "Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\"")
				.build();
	}

	@Override
	public Map<String, String> getFormatFields() {
		return new MapBuilder<String, String>()
				.put("DP", "Number=1,Type=Integer,Description=\"# high-quality bases\"")
				.put("DV", "Number=1,Type=Integer,Description=\"# high-quality non-reference bases\"")
				.build();
	}
}
