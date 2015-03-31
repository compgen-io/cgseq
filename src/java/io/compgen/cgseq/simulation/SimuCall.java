package io.compgen.cgseq.simulation;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.common.IterUtils;
import io.compgen.common.IterUtils.MapFunc;
import io.compgen.common.StringUtils;
import io.compgen.ngsutils.pileup.PileupRecord;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCall;

import java.util.List;
import java.util.Random;

/**
 * Simulates base calls for testing variant callers
 * 
 * This is effectively simulates a bionomial distribution, allowing for a
 * configurable variant allele frequency (in bionomial freq = 0.5). It also
 * allows for random sequencing error 
 * 
 * @author mbreese
 *
 */

@Command(name = "simucall")
public class SimuCall {
	private double variantAlleleFreq = 0.5;
	private double strandFreq = 0.5;
	private int phredScore = 30;
	private int depth = 30;
	
	private int times = 1;
	
	private String majorCall = "A";
	private String altCall = "C";
	private String errorCall = "G";

	@Exec
	public void exec() {
		for (int i=0; i<times; i++) {
			List<PileupBaseCall> calls = simulate();
			System.out.println(StringUtils.join("", IterUtils.map(calls,  new MapFunc<PileupBaseCall, String>(){
				@Override
				public String map(PileupBaseCall call) {
					if (call.plusStrand) {
						return call.call.toUpperCase();
					}
					return call.call.toLowerCase();
				}})));
		}
	}

	public List<PileupBaseCall> simulate() {
		double errorRate = Math.pow(10, phredScore / -10);
//		System.err.println("Error rate = "+ errorRate + " (phred="+phredScore+")");
		
		Random rand = new Random();
		String calls = "";
		String quals = "";
		
		for (int i=0; i<depth; i++) {
			quals += (char)(phredScore + 33);

			boolean plusStrand = (rand.nextDouble() > strandFreq);
			
			if (rand.nextDouble() < errorRate) {
				if (plusStrand) {
					calls += errorCall.toUpperCase();
				} else {
					calls += errorCall.toLowerCase();
				}
			} else if (rand.nextDouble() < variantAlleleFreq) {
				// yes, another call to nextDouble() -> we need to split
				// the error vs. alt-call test
				if (plusStrand) {
					calls += altCall.toUpperCase();
				} else {
					calls += altCall.toLowerCase();
				}
			} else {
				if (plusStrand) {
					calls += ".";
				} else {
					calls += ",";
				}
			}
		}
		
		String pileupLine = "chrTest\t101\t"+majorCall+"\t"+depth+"\t"+calls+"\t"+quals+"\n";
		return PileupRecord.parse(pileupLine).getSampleRecords(0).calls;

	}

	@UnnamedArg(name="depth")
	public void setDepth(int depth) {
		this.depth = depth;
	}

	@Option(name="strand-freq")
	public void setStrandFreq(double strandFreq) {
		this.strandFreq = strandFreq;
	}

	@Option(name="alt-freq")
	public void setVariantAlleleFreq(double variantAlleleFreq) {
		this.variantAlleleFreq = variantAlleleFreq;
	}

	@Option(name="times")
	public void setTimes(int times) {
		this.times = times;
	}

	@Option(name="phred")
	public void setPhredScore(int phredScore) {
		this.phredScore = phredScore;
	}

	@Option(name="major")
	public void setMajorCall(String majorCall) {
		this.majorCall = majorCall;
	}

	@Option(name="alt")
	public void setAltCall(String altCall) {
		this.altCall = altCall;
	}

	@Option(name="error")
	public void setErrorCall(String errorCall) {
		this.errorCall = errorCall;
	}
	
}
