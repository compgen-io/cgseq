package io.compgen.cgseq.cli.copynumber;

import io.compgen.cgseq.CGSeq;
import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.ListBuilder;
import io.compgen.common.StringUtils;
import io.compgen.common.TabWriter;
import io.compgen.ngsutils.pileup.PileupReader;
import io.compgen.ngsutils.pileup.PileupRecord;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.stat.regression.SimpleRegression;

@Command(name="breakpoint", desc="Find likely breakpoints across a genome", category="copy-number")
public class BreakpointFinder extends AbstractOutputCommand {
	
	public class PileupCount{
		public final String ref;
		public final int pos;
		public final int count;
		
		public PileupCount(PileupRecord record) {
			this.ref = record.ref;
			this.pos = record.pos;
			this.count = record.getSampleCount(0);			
		}
		public PileupCount(String ref, int pos, int count) {
			this.ref = ref;
			this.pos = pos;
			this.count = count;
		}
	}
	
	private String filename = "-";
	private int windowSize = 10000;
	private int stepSize = 2000;
	private double fractionIgnore = 0.1;
	private double fractionRegression = 0.2;

	private int minBaseQual = 13;
	private int minMapQ = 0;
	private boolean properPairs = false;

    @Option(desc="Only count properly-paired reads", name="paired")
    public void setProperPairs(boolean properPairs) {
    	this.properPairs = properPairs;
    }

    @Option(desc="Minimum base quality", name="min-basequal", defaultValue="13")
    public void setMinBaseQual(int minBaseQual) {
    	this.minBaseQual = minBaseQual;
    }

    @Option(desc="Minimum read mapping quality (MAPQ)", name="min-mapq", defaultValue="0")
    public void setMinMapQual(int minMapQ) {
    	this.minMapQ = minMapQ;
    }
	
    @Option(desc="Fraction of counts to ignore (0.0-1.0)", name="fraction-ignore", defaultValue="0.1")
    public void setFractionIgnore(double fractionIgnore) {
    	this.fractionIgnore = fractionIgnore;
    }

    @Option(desc="Fraction of counts to use in regression (0.0-1.0)", name="fraction-regression", defaultValue="0.2")
    public void setFractionRegression(double fractionRegression) {
    	this.fractionRegression = fractionRegression;
    }

    @Option(desc="Window size (bp)", name="window-size", defaultValue="10000")
    public void setWindowSize(int windowSize) {
    	this.windowSize = windowSize;
    }

    @Option(desc="Step size (bp)", name="step-size", defaultValue="2000")
    public void setStepSize(int stepSize) {
    	this.stepSize = stepSize;
    }

    @UnnamedArg(name = "BAM")
    public void setFilename(String filename) {
        this.filename = filename;
    }

	public BreakpointFinder() {
	}


	@Exec
	public void exec() throws Exception {
		TabWriter writer = new TabWriter(out);
        writer.write_line("## program: " + CGSeq.getVersion());
        writer.write_line("## cmd: " + CGSeq.getArgs());
        writer.write_line("## filename: " + filename);
        writer.write_line("## min-mapq: " + minMapQ);
        writer.write_line("## min-base-qual: " + minBaseQual);
        writer.write_line("## proper-pairs: " + properPairs);
		writer.write("chrom", "start", "end", "median", "slope", "intercept", "mse", "rsquared", "mse_median");
		writer.eol();

		ListBuilder<String> lb = new ListBuilder<String>("samtools", "mpileup");
		if (minMapQ > -1) {
			lb.add("-q", ""+minMapQ);
		}
		if (minBaseQual > -1) {
			lb.add("-Q", ""+minBaseQual);
		}
		if (properPairs) {
			lb.add("--rf", "2");
		}
		
		lb.add(filename);
		
		ProcessBuilder pb = new ProcessBuilder(lb.list());
		if (verbose) {
			System.err.println(StringUtils.join(" ", pb.command()));
		}
		Process proc = pb.start();
		InputStream bis = new BufferedInputStream(proc.getInputStream());
		PileupReader reader = new PileupReader(bis);

		String currentChrom = null;
		int curStart = 0;

		List<PileupCount> buffer = new ArrayList<PileupCount>();
		Iterator<PileupRecord> it = reader.iterator();

		PileupRecord record = null;
		while (it.hasNext()) {
			if (record == null) {
				 record = it.next();
			}

			if (currentChrom == null || !record.ref.equals(currentChrom)) {
				if (currentChrom!=null) {
					calculateRegression(buffer, writer, curStart, curStart + windowSize);
					buffer.clear();
				}
				currentChrom = record.ref;
				curStart = 0;

			} else if (record.pos > curStart + windowSize) {
				if (buffer.size() > 0) {
					calculateRegression(buffer, writer, curStart, curStart + windowSize);
				}

				curStart = curStart + stepSize;
				
				while (buffer.size() > 0 && buffer.get(0).pos < curStart) {
					buffer.remove(0);
				}
				
			} else {
				buffer.add(new PileupCount(record));
				record = null;
			}
		}

		if (buffer.size() > 0) {
			calculateRegression(buffer, writer, curStart, curStart + windowSize);
		}
		
		writer.close();
		reader.close();
	}
	
	private void calculateRegression(List<PileupCount> buffer, TabWriter writer, int start, int end) throws IOException {
		int[] counts = new int[buffer.size()];
		for (int i=0; i<counts.length; i++) {
			counts[i] = buffer.get(i).count;
		}
		
        int skip = (int) (counts.length*fractionIgnore);

//      int low = skip;
//		int high = counts.length - skip;
//		int quarterLow = ((mid-low) / 2) + low;
//		int quarterhigh = ((high-mid) / 2) + mid;
		
		Arrays.sort(counts);
		
		// first quartile
		
		int mid = counts.length / 2;
		int split = (int) (counts.length * fractionRegression * 0.5);
		
        SimpleRegression first = calcRegression(counts, mid-split, mid+split);

		double mse = calculateMSE(first, counts, skip);
		double rsquare = first.getRSquare();
		double slope = first.getSlope();
		double intercept = first.getIntercept();
		int median = counts[(int) (counts.length*0.5)];
		
//		int quarter = 1;
//
//        double mse = calculateMSE(second, counts, skip);
//        if (mse > worstMSE) {
//    		worstMSE = mse;
//    		rsquare = second.getRSquare();
//    		slope = second.getSlope();
//    		intercept = second.getIntercept();
//    		quarter = 2;
//        }
//
//        mse = calculateMSE(third, counts, skip);
//        if (mse > worstMSE) {
//    		worstMSE = mse;
//    		rsquare = third.getRSquare();
//    		slope = third.getSlope();
//    		intercept = third.getIntercept();
//    		quarter = 3;
//        }
//
//        mse = calculateMSE(fourth, counts, skip);
//        if (mse > worstMSE) {
//    		worstMSE = mse;
//    		rsquare = fourth.getRSquare();
//    		slope = fourth.getSlope();
//    		intercept = fourth.getIntercept();
//    		quarter = 4;
//        }
        
		writer.write(buffer.get(0).ref);
		writer.write(start);
		writer.write(end);
		writer.write(median);
		writer.write(slope);
		writer.write(intercept);
		writer.write(mse);
		writer.write(rsquare);
		writer.write(mse/median);
		writer.eol();

        
	}
	
	private SimpleRegression calcRegression(int[] counts, int low, int high) {
		SimpleRegression reg = new SimpleRegression();
        for (int i=low; i<high; i++) {
            reg.addData(i,counts[i]);
        }
        return reg;
	}
	
	private double calculateMSE(SimpleRegression reg, int[] counts, int skip) {
		double acc = 0.0;
		
		for (int i=skip; i<counts.length-skip; i++) {
			acc += Math.pow(reg.predict(i) - counts[i],2);
		}
		
		return acc / counts.length;
	}
	
//
//	private double calcCummulativeDistance(int[] ar) {
//		int total = 0;
//		for (int i=0; i<ar.length;i++) {
//			total += ar[i];
//		}
//
//		double step = 1.0 / ar.length;
//		double idealAcc = 0.0;
//		double sampleAcc = 0.0;
//		
//		double dist = 0.0;
//		for (int i=0; i<ar.length;i++) {
//			sampleAcc += ar[i];
//			idealAcc += step;
//			dist += Math.abs((sampleAcc/total) - idealAcc);
//		}
//
//		return dist;
//	}
//
//	private double calcTumorNormalDiff(int[] normal, int[] tumor) {
//		int normalTotal = 0;
//		int tumorTotal = 0;
//
//		for (int i=0; i<normal.length;i++) {
//			normalTotal += normal[i];
//			tumorTotal += tumor[i];
//		}
//
//		double normalAcc = 0.0;
//		double tumorAcc = 0.0;
//		
//		double diff = 0.0;
//
//		for (int i=0; i<normal.length;i++) {
//			normalAcc += normal[i];
//			tumorAcc += tumor[i];
//			
//			diff += Math.abs((normalAcc/normalTotal) - (tumorAcc/tumorTotal));
//		}
//
//		return diff / normal.length;
//	}

}
