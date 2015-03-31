package org.ngsutils.cgsutils.cli.copynumber;

import java.io.BufferedInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.ngsutils.annotation.GenomeRegion;
import org.ngsutils.cgsutils.CGSUtils;
import org.ngsutils.cgsutils.support.pileup.PileupReader;
import org.ngsutils.cgsutils.support.pileup.PileupRecord;
import org.ngsutils.cli.AbstractOutputCommand;
import org.ngsutils.cli.Command;
import org.ngsutils.support.ListUtils;
import org.ngsutils.support.StringLineReader;
import org.ngsutils.support.StringUtils;
import org.ngsutils.support.TabWriter;
import org.ngsutils.support.stats.StatUtils;

import com.lexicalscope.jewel.cli.ArgumentValidationException;
import com.lexicalscope.jewel.cli.CommandLineInterface;
import com.lexicalscope.jewel.cli.Option;

@CommandLineInterface(application="cgsutils copynumber")
@Command(name="copynumber", 
		 desc="Find the overall copy-number for the tumor sample across the a region (mpileup input, NT).", 
		 cat="copy-number", 
		 doc="You can use either mpileup as an input or BAM files. If you use mpileup (Normal/Tumor order),\n"
		 		+ "then the copy-number will be calculated across the whole mpileup region. Otherwise, you\n"
		 		+ "can also use normal and tumor BAM files with a BED file to define the regions used to\n"
		 		+ "calculate copy numbers. In this case, 'samtools' must be present in the $PATH.\n\n"
		 		+ "For each region, the output will be a tab-delimited line with the tumor/normal ratio (log2)\n"
		 		+ "and an estimated copy-number. If --norm-total and --tumor-total are included, these will be\n"
		 		+ "used to normalize the counts to give a more accurate copy number estimate. The greater the\n"
		 		+ "size of the region, the more accurate the estimate will be (assuming that there aren't any\n"
		 		+ "copy-number changes within the region).\n\n"
		 		+ "The log-ratio is determined by taking the log2 ratio of tumor/normal using the the median\n"
		 		+ "read count. Copy number is then 2^(log2-ratio + 1)."
		 )

public class PileupCopyNumber extends AbstractOutputCommand {
	public class CopyNumberRecord {
		public final String chrom;
		public final int start;
		public final int end;
		public final double ratio;
		public final double copyNumber;

		public CopyNumberRecord(String chrom, int start, int end, double ratio, double copyNumber) {
			this.chrom = chrom;
			this.start = start;
			this.end = end;
			this.ratio = ratio;
			this.copyNumber = copyNumber;
		}
	}

	private String pileupFilename = null;
	private String normalFilename = null;
	private String tumorFilename = null;
	private String bedFilename = null;
	private String region = null;
	private int normalTotal = -1;
	private int tumorTotal = -1;

    @Option(description="Normal BAM file", longName="norm", defaultToNull=true)
    public void setNormalFilename(String filename) {
    	this.normalFilename = filename;
    }

    @Option(description="Tumor BAM file", longName="tumor", defaultToNull=true)
    public void setTumorFilename(String filename) {
    	this.tumorFilename = filename;
    }
    
    @Option(description="Regions BED file", longName="bed", defaultToNull=true)
    public void setBEDFilename(String filename) {
    	this.bedFilename = filename;
    }
    @Option(description="Region to find copy-number for (using BAM files)", longName="region", defaultToNull=true)
    public void setRegion(String region) {
    	this.region = region;
    }
    
    @Option(description="Normal total read count (optional)", longName="norm-total", defaultValue="-1")
    public void setNormalTotal(int count) {
    	this.normalTotal = count;
    }

    @Option(description="Tumor total read count (optional)", longName="tumor-total", defaultValue="-1")
    public void setTumorTotal(int count) {
        	this.tumorTotal = count;
    }

    @Option(description="Pileup file", longName="pileup", defaultToNull=true)
    public void setPileupFilename(String pileupFilename) {
        this.pileupFilename = pileupFilename;
    }

	public PileupCopyNumber() {
	}

	public static final double log2(double d) {
		return Math.log(d) / Math.log(2);
	}
	
	@Override
	public void exec() throws Exception {
		if (pileupFilename == null && (tumorFilename == null || normalFilename == null || (region == null && bedFilename == null))) {
			throw new ArgumentValidationException("You must specify either an mpileup file (or stdin), or a normal BAM file, a tumor BAM file, and a region/BED file!");
		}
		
		TabWriter writer = new TabWriter(out);
        writer.write_line("## program: " + CGSUtils.getVersion());
        writer.write_line("## cmd: " + CGSUtils.getArgs());
		if (pileupFilename != null) {
			writer.write_line("## pileup-input: " + pileupFilename);
		} else {
			writer.write_line("## bed-regions: " + bedFilename);
			writer.write_line("## normal: " + normalFilename);
			writer.write_line("## tumor: " + tumorFilename);
		}
        writer.write_line("## tumor-total: " + tumorTotal);
        writer.write_line("## normal-total: " + normalTotal);
		
		if (pileupFilename != null) {
			writer.write("chrom", "start", "end", "ratio (log2)", "copy-number");
			writer.eol();
			PileupReader reader = new PileupReader(pileupFilename);
			CopyNumberRecord record = calcCopyNumber(reader, normalTotal, tumorTotal);
			writer.write(record.chrom);
			writer.write(record.start);
			writer.write(record.end);
			writer.write(record.ratio);
			writer.write(record.copyNumber);
			writer.eol();
			reader.close();
		} else if (region == null){
			writer.write("name", "chrom", "start", "end", "ratio (log2)", "copy-number");
			writer.eol();
			
			StringLineReader strReader = new StringLineReader(bedFilename);
			for (String line: strReader) {
				String[] cols = StringUtils.strip(line).split("\t");
				String chrom = cols[0];
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				String name = cols[3];

				if (verbose) {
					System.err.println(name);
				}
				
				ProcessBuilder pb = new ProcessBuilder("samtools", "mpileup", "-r", chrom + ":" + ""+(start+1) + "-" + ""+end , normalFilename,tumorFilename);
				if (verbose) {
					System.err.println(StringUtils.join(" ", pb.command()));
				}
				Process proc = pb.start();
				InputStream bis = new BufferedInputStream(proc.getInputStream());
				PileupReader reader = new PileupReader(bis);

				CopyNumberRecord record = calcCopyNumber(reader, normalTotal, tumorTotal);
				writer.write(name);
				writer.write(chrom);
				writer.write(start);
				writer.write(end);
				writer.write(record.ratio);
				writer.write(record.copyNumber);
				writer.eol();
				reader.close();
				proc.waitFor();
				proc.getErrorStream().close();
			    proc.getInputStream().close();
			    proc.getOutputStream().close();
			    proc.destroy();

			}
			strReader.close();
		} else {
			writer.write("chrom", "start", "end", "ratio (log2)", "copy-number");
			writer.eol();
			
			ProcessBuilder pb = new ProcessBuilder("samtools", "mpileup", "-r", region , normalFilename,tumorFilename);
			if (verbose) {
				System.err.println(StringUtils.join(" ", pb.command()));
			}
			
			Process proc = pb.start();
			InputStream bis = new BufferedInputStream(proc.getInputStream());
			PileupReader reader = new PileupReader(bis);

			GenomeRegion gen = GenomeRegion.parse(region, true);
			
			CopyNumberRecord record = calcCopyNumber(reader, normalTotal, tumorTotal);
			writer.write(gen.ref);
			writer.write(gen.start);
			writer.write(gen.end);
			writer.write(record.ratio);
			writer.write(record.copyNumber);
			writer.eol();
			reader.close();
			proc.waitFor();
			proc.getErrorStream().close();
		    proc.getInputStream().close();
		    proc.getOutputStream().close();
		    proc.destroy();
		}
		writer.close();
	}

	public CopyNumberRecord calcCopyNumber(PileupReader reader, int normalTotal, int tumorTotal) {
		List<Integer> normalCounts = new ArrayList<Integer>(); 
		List<Integer> tumorCounts = new ArrayList<Integer>(); 
		
		String chrom = null;
		int start = -1;
		int end = -1;

		for (PileupRecord pileup: reader) {
			
			if (chrom == null) {
				chrom = pileup.ref;
				start = pileup.pos;
			}
			
			end = pileup.pos;
			
			normalCounts.add(pileup.getSampleCount(0));
			tumorCounts.add(pileup.getSampleCount(1));
		}

		int[] norm = ListUtils.intListToArray(normalCounts);
		int[] tumor = ListUtils.intListToArray(tumorCounts);

		double medianRatio = calcCopyRatio(norm, tumor, normalTotal, tumorTotal);
		
		double copyNumber = calcCopyNumber(medianRatio);
		
		return new CopyNumberRecord(chrom, start, end, medianRatio, copyNumber);

	}
	
	public static double calcCopyNumber(double ratio) {
		return Math.pow(2, ratio + 1);
	}
	public static double calcCopyRatio(int[] norm, int[] tumor) {
		return calcCopyRatio(norm, tumor, -1, -1);
	}
	
	public static double calcCopyRatio(int[] norm, int[] tumor, int normNormalization, int tumorNormalization) {
		double[] normalizedNorm;
		double[] normalizedTumor;

		if (normNormalization > 0 && tumorNormalization > 0) {
			normalizedNorm = normalizeLog2(norm, normNormalization);
			normalizedTumor = normalizeLog2(tumor, tumorNormalization);
		} else {
			normalizedNorm = StatUtils.log2(norm);
			normalizedTumor = StatUtils.log2(tumor);
		}
		
		Arrays.sort(normalizedNorm);
		Arrays.sort(normalizedTumor);
		
		double normMedian = StatUtils.median(normalizedNorm);
		double tumorMedian = StatUtils.median(normalizedTumor);
		
		return tumorMedian - normMedian;
		
//		double[] ratio = ListUtils.doubleArrayPairMap(normalizedNorm, normalizedTumor, new PairExec<Double>() {
//			@Override
//			public Double func(Double one, Double two) {
//				// Tumor - Normal is a log-ratio
//				return two - one;
//			}});
//
//		Arrays.sort(ratio);
//		return StatUtils.median(ratio);
	}
	 
	public static double[] normalizeLog2(int[] vals, int total) {
		double[] out = new double[vals.length];
		double totalLog2 = StatUtils.log2(total);
		for (int i=0; i<vals.length; i++) {
			out[i] = StatUtils.log2(vals[i]) - totalLog2;
		}
		return out;
	}
	
}
