package io.compgen.cgseq.cli.copynumber;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import io.compgen.cgseq.CGSeq;
import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.IterUtils;
import io.compgen.common.StringUtils;
import io.compgen.common.TabWriter;
import io.compgen.common.progress.ProgressMessage;
import io.compgen.common.progress.ProgressStats;
import io.compgen.common.progress.ProgressUtils;
import io.compgen.ngsutils.annotation.GenomeSpan;
import io.compgen.ngsutils.pileup.BAMPileup;
import io.compgen.ngsutils.pileup.PileupRecord;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@Command(name="bp-dist", desc="Calculate somatic/germline distance across sliding windows", category="copy-number")
public class BreakpointFinder extends AbstractOutputCommand {
	
	public class BPWindowStats{
		public final double dist;
		public final int maxPos;
		
		public BPWindowStats(double dist, int maxPos) {
			this.dist = dist;
			this.maxPos = maxPos;
		}
	}
	
	private String germlineFname=null;
	private String somaticFname=null;
	
	private int windowSize = 10000;
	private int stepSize = 2500;

	private int minBaseQual = 13;
	private int minMapQ = 0;
	private boolean properPairs = false;

	private int filterFlags = 0;
    private int requiredFlags = 0;

    private String region = null;
    
    @Option(desc="Only calculated breakpoints for this region", name="region", charName="R")
    public void setRegion(String region) {
    	this.region = region;
    }

    @Option(desc="Only count properly-paired reads", name="paired")
    public void setProperPairs(boolean properPairs) {
    	this.properPairs = properPairs;
    	this.requiredFlags |= 0x2;
    }

    @Option(desc="Minimum base-quality score", name="min-basequal", defaultValue="30")
    public void setMinBaseQual(int minBaseQual) {
    	this.minBaseQual = minBaseQual;
    }

    @Option(desc="Minimum alignment mapping score (MAPQ)", name="min-mapq", defaultValue="10")
    public void setMinMapQual(int minMapQ) {
    	this.minMapQ = minMapQ;
    }

    @Option(desc="Filter flags", name="filter-flags")
    public void setFilterFlags(int filterFlags) {
    	this.filterFlags = filterFlags;
    }
    
    @Option(desc="Required flags", name="required-flags")
    public void setRequiredFlags(int requiredFlags) {
    	this.requiredFlags = requiredFlags;
    }
    
    @Option(desc="Window size (bp)", name="window-size", defaultValue="10000")
    public void setWindowSize(int windowSize) {
    	this.windowSize = windowSize;
    }

    @Option(desc="Step size (bp)", name="step-size", defaultValue="5000")
    public void setStepSize(int stepSize) {
    	this.stepSize = stepSize;
    }

    @UnnamedArg(name = "GERMLINE SOMATIC")
    public void setFilename(String[] filenames) throws CommandArgumentException {
        if (filenames.length!=2) {
        	throw new CommandArgumentException("You must specify both a germline and somatic sample");
        }
        
        germlineFname = filenames[0];
        somaticFname = filenames[1];
        
    }

	public BreakpointFinder() {
	}


	@Exec
	public void exec() throws CommandArgumentException, IOException {
		BAMPileup pileup = new BAMPileup(germlineFname, somaticFname);
		pileup.setDisableBAQ(true);
		pileup.setExtendedBAQ(false);
		pileup.setFlagFilter(filterFlags);
		pileup.setFlagRequired(requiredFlags);
		pileup.setMinBaseQual(minBaseQual);
		pileup.setMinMappingQual(minMapQ);

		SamReader bam = SamReaderFactory.makeDefault().open(new File(germlineFname));
		final SAMFileHeader header = bam.getFileHeader();

		GenomeSpan regionSpan = null;
		if (region != null) {
			if (region.indexOf(':') > -1) {
				regionSpan = GenomeSpan.parse(region);
				if (header.getSequence(regionSpan.ref) == null) {
					throw new CommandArgumentException("Region: "+ region+" not found in this BAM file!");
				}

			} else {
				// this is just a raw chrom, we need to pull the length 
				if (header.getSequence(region) == null) {
					throw new CommandArgumentException("Region: "+ region+" not found in this BAM file!");
				}
				region = region+":1-"+header.getSequence(region).getSequenceLength();
				regionSpan = GenomeSpan.parse(region);
			}

		}

		final long totalGenomeSize;

		if (regionSpan == null) {
			long tmp = 0;
			for (SAMSequenceRecord seq: header.getSequenceDictionary().getSequences()) {
				tmp += seq.getSequenceLength();
			}
			totalGenomeSize = tmp;
		} else {
			totalGenomeSize = regionSpan.size();
		}

		
		TabWriter writer = new TabWriter(out);
        writer.write_line("## program: " + CGSeq.getVersion());
        writer.write_line("## cmd: " + CGSeq.getArgs());
        writer.write_line("## germline: " + germlineFname);
        writer.write_line("## somatic: " + somaticFname);
        writer.write_line("## min-mapq: " + minMapQ);
        writer.write_line("## min-base-qual: " + minBaseQual);
        writer.write_line("## proper-pairs: " + properPairs);
		writer.write_line("## pileupCommand="+StringUtils.join(" ", pileup.getCommand(regionSpan)));

		for (SAMSequenceRecord seq: header.getSequenceDictionary().getSequences()) {
	        writer.write_line("## ref "+seq.getSequenceName()+" " + seq.getSequenceLength());
		}

		
		String currentChrom = null;
		int curStart = 0;

		
		final long[] progressPos = new long[] {0l,0l}; 
 		
		final List<PileupRecord> buffer = new ArrayList<PileupRecord>();
//		final Map<BPPos,Double> stats = new HashMap<BPPos, Double>();

		
		Iterator<PileupRecord> it = ProgressUtils.getIterator(new File(germlineFname).getName()+" / " + new File(somaticFname).getName(), pileup.pileup(regionSpan), new ProgressStats(){
			@Override
			public long size() {
				return totalGenomeSize;
			}
			@Override
			public long position() {
				return progressPos[0] + progressPos[1];
			}}, 
			
			new ProgressMessage<PileupRecord>(){
				@Override
				public String msg(PileupRecord current) {
						return current.ref+":"+current.pos;
				}});
		
		
		for (PileupRecord record: IterUtils.wrap(it)) {
			if (currentChrom == null || !record.ref.equals(currentChrom)) {
				if (buffer.size() > 0) {
					BPWindowStats diff = calcCumulativeDistance(buffer);
//					stats.put(diff.pos, diff.dist);
					writer.write(currentChrom);
					writer.write(curStart);
					writer.write(buffer.get(buffer.size()-1).pos);
					writer.write(Double.isNaN(diff.dist) ? "":""+diff.dist);
					writer.write(diff.maxPos < 0 ? "": ""+diff.maxPos);
					writer.eol();
					buffer.clear();
					
					progressPos[0] += header.getSequence(currentChrom).getSequenceLength();
					progressPos[1] = 0;
					
				}
				currentChrom = record.ref;
				curStart = 0;
			} 
			if (record.pos > curStart + windowSize) {
				if (buffer.size() > 0) {
					BPWindowStats diff = calcCumulativeDistance(buffer);
//					stats.put(diff.pos, diff.dist);
					writer.write(currentChrom);
					writer.write(curStart);
					writer.write(curStart + windowSize);
					writer.write(Double.isNaN(diff.dist) ? "":""+diff.dist);
					writer.write(diff.maxPos < 0 ? "": ""+diff.maxPos);
					writer.eol();
				}
				while (record.pos > curStart + windowSize) {
					curStart += stepSize;
				}
				
				while (buffer.size() > 0 && buffer.get(0).pos < curStart) {
					buffer.remove(0);
				}
			}

			buffer.add(record);
			progressPos[1] = record.pos;
		}

		if (buffer.size() > 0) {
			BPWindowStats diff = calcCumulativeDistance(buffer);
//			stats.put(diff.pos, diff.dist);
			writer.write(currentChrom);
			writer.write(curStart);
			writer.write(buffer.get(buffer.size()-1).pos);
			writer.write(Double.isNaN(diff.dist) ? "":""+diff.dist);
			writer.write(diff.maxPos < 0 ? "": ""+diff.maxPos);
			writer.eol();
			buffer.clear();
		}
		
//		System.err.println("Calculating mean/stddev of all distances");
//		double[] distances = new double[stats.size()];
//		int i=0;
//		for (Double d: stats.values()) {
//			distances[i++] = d;
//		}
//		MeanStdDev msd = StatUtils.calcMeanStdDev(distances);
//		System.err.println("Mean  : " + msd.mean);
//		System.err.println("StdDev: " + msd.stddev);
//		
//		double thres = (msd.stddev*2) + msd.mean;
//		System.err.println("Threshold (2-sigma): " + thres);
//
//		writer.write_line("#window-count: " + distances.length);
//		writer.write_line("#mean: " + msd.mean);
//		writer.write_line("#stddev: " + msd.stddev);
//		writer.write_line("#threshold (2-sigma): " + thres);
//
//		writer.write("chrom", "pos", "cumulative_dist");
//		writer.eol();
//
//		for (Entry<BPPos, Double> ev: stats.entrySet()) {
//			if (ev.getValue() > thres) {
//				writer.write(ev.getKey().chrom);
//				writer.write(ev.getKey().pos);
//				writer.eol();
//			}
//		}
		
		writer.close();
	}
	
	private BPWindowStats calcCumulativeDistance(List<PileupRecord> buffer) {
		long germTotal = 0;
		long somTotal = 0;

		for (int i=0; i<buffer.size(); i++) {
			germTotal += buffer.get(i).getSampleCount(0);
			somTotal += buffer.get(i).getSampleCount(1);
		}

		double germAcc = 0.0;
		double somAcc = 0.0;
		
		double diffAcc = 0.0;
		double maxDiff = -1.0;
		int maxPos = -1;

		for (int i=0; i<buffer.size(); i++) {
			germAcc = buffer.get(i).getSampleCount(0);
			somAcc = buffer.get(i).getSampleCount(1);
			
			double diff = Math.abs((germAcc/germTotal) - (somAcc/somTotal));
			
			if (diff > maxDiff) {
				maxPos = buffer.get(i).pos;
				maxDiff = diff;
			}
			
			diffAcc += diff;

		}

		return new BPWindowStats(diffAcc / buffer.size(), maxPos);
	}

}
