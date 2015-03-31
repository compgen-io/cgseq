package io.compgen.cgseq.cli.copynumber;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.TabWriter;
import io.compgen.ngsutils.pileup.PileupReader;
import io.compgen.ngsutils.pileup.SlidingWindowPileup;

@Command(name="breakpoints", desc="Find likely breakpoints across a chromosome (mpileup input, NT)", category="copy-number")
public class Breakpoints extends AbstractOutputCommand {
	private String filename = "-";
	private int windowSize = 10000;
	private int stepSize = 2000;

    @Option(desc="Window size (default: 10000)", name="window-size", defaultValue="10000")
    public void setWindowSize(int windowSize) {
    	this.windowSize = windowSize;
    }

    @Option(desc="Step size (default: 2000)", name="step-size", defaultValue="2000")
    public void setStepSize(int stepSize) {
    	this.stepSize = stepSize;
    }

    @UnnamedArg(name = "FILE")
    public void setFilename(String filename) {
        this.filename = filename;
    }

	public Breakpoints() {
	}

	public static final double log2(double d) {
		return Math.log(d) / Math.log(2);
	}
	
	@Exec
	public void exec() throws Exception {
		TabWriter writer = new TabWriter(out);
		writer.write("chrom", "start", "end", "normal_dist", "tumor_dist", "ratio (log2)", "tumor_normal_diff");
		writer.eol();

		PileupReader reader = new PileupReader(filename);
		
		int i=0;
		
		for (SlidingWindowPileup counts:SlidingWindowPileup.readMPileup(reader, windowSize, stepSize)) {
			if (verbose && (i++ % 100 == 0)) {
				System.err.println(counts.ref+":"+counts.start);
			}
			writer.write(counts.ref, ""+counts.start, ""+counts.end);
			double normalDist = calcCummulativeDistance(counts.normalCounts);
			double tumorDist = calcCummulativeDistance(counts.tumorCounts);
			writer.write(normalDist);
			writer.write(tumorDist);
			writer.write(log2(tumorDist) - log2(normalDist));
			writer.write(calcTumorNormalDiff(counts.normalCounts, counts.tumorCounts));
			writer.eol();
		}
		writer.close();
		reader.close();
	}
	
	private double calcCummulativeDistance(int[] ar) {
		int total = 0;
		for (int i=0; i<ar.length;i++) {
			total += ar[i];
		}

		double step = 1.0 / ar.length;
		double idealAcc = 0.0;
		double sampleAcc = 0.0;
		
		double dist = 0.0;
		for (int i=0; i<ar.length;i++) {
			sampleAcc += ar[i];
			idealAcc += step;
			dist += Math.abs((sampleAcc/total) - idealAcc);
		}

		return dist;
	}

	private double calcTumorNormalDiff(int[] normal, int[] tumor) {
		int normalTotal = 0;
		int tumorTotal = 0;

		for (int i=0; i<normal.length;i++) {
			normalTotal += normal[i];
			tumorTotal += tumor[i];
		}

		double normalAcc = 0.0;
		double tumorAcc = 0.0;
		
		double diff = 0.0;

		for (int i=0; i<normal.length;i++) {
			normalAcc += normal[i];
			tumorAcc += tumor[i];
			
			diff += Math.abs((normalAcc/normalTotal) - (tumorAcc/tumorTotal));
		}

		return diff / normal.length;
	}

}
