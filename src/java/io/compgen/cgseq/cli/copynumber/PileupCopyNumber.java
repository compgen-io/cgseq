package io.compgen.cgseq.cli.copynumber;

import io.compgen.cgseq.CGSeq;
import io.compgen.cgseq.support.MapCount;
import io.compgen.cgseq.variant.SkellamVariantCaller;
import io.compgen.cgseq.variant.VariantResults;
import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.ListBuilder;
import io.compgen.common.Pair;
import io.compgen.common.StringLineReader;
import io.compgen.common.StringUtils;
import io.compgen.common.TabWriter;
import io.compgen.ngsutils.annotation.GenomeSpan;
import io.compgen.ngsutils.pileup.PileupReader;
import io.compgen.ngsutils.pileup.PileupRecord;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCall;
import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCallOp;
import io.compgen.ngsutils.pileup.PileupRecord.PileupSampleRecord;
import io.compgen.ngsutils.support.stats.StatUtils;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Command(name="copynumber", 
		 desc="Find the overall copy-number for the tumor sample across the a region (mpileup input, NT).", 
		 category="copy-number", 
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
		 		+ "read count. Copy number is then 2^(log2-ratio + 1).\n\n"
		 		+ "Average minor-allele frequency is calculated based on the minor-allele frequency for *only*\n"
		 		+ "heterozygous germline-variants. Germline variants calculated using the Skellam distribution."
		 )

public class PileupCopyNumber extends AbstractOutputCommand {
	public class CopyNumberRecord {
		public final String chrom;
		public final int start;
		public final int end;
		public final double ratio;
		public final double copyNumber;
		public final double normMedian;
		public final double tumorMedian;
		public final double mafAve;
		public final double mafStdev;
		public final int mafCount;
		

		public CopyNumberRecord(String chrom, int start, int end, double ratio, double copyNumber, double normMedian, double tumorMedian, double mafAve, double mafStdev, int mafCount) {
			this.chrom = chrom;
			this.start = start;
			this.end = end;
			this.ratio = ratio;
			this.copyNumber = copyNumber;
			this.normMedian = normMedian;
			this.tumorMedian = tumorMedian;
			this.mafAve = mafAve;
			this.mafStdev = mafStdev;
			this.mafCount = mafCount;
		}
	}

	private String pileupFilename = null;
	private String normalFilename = null;
	private String tumorFilename = null;
	private String bedFilename = null;
	private String region = null;
	private int normalTotal = -1;
	private int tumorTotal = -1;

	private int minMAFDepth = 10;
	private int minBaseQual = 13;
	private int minMapQ = 0;

	private boolean properPairs = false;
	private SkellamVariantCaller germlineCaller = null;

    @Option(desc="Only count properly-paired reads", name="paired")
    public void setProperPairs(boolean properPairs) {
    	this.properPairs = properPairs;
    }
    
    @Option(desc="Minimum depth for MAF calculation", name="min-maf-depth", defaultValue="10")
    public void setMinMAFDepth(int minMAFDepth) {
    	this.minMAFDepth = minMAFDepth;
    }

    @Option(desc="Minimum base quality", name="min-basequal", defaultValue="13")
    public void setMinBaseQual(int minBaseQual) {
    	this.minBaseQual = minBaseQual;
    }

    @Option(desc="Minimum read mapping quality (MAPQ)", name="min-mapq", defaultValue="0")
    public void setMinMapQual(int minMapQ) {
    	this.minMapQ = minMapQ;
    }
	
    @Option(desc="Normal BAM file", name="norm", helpValue="fname")
    public void setNormalFilename(String filename) {
    	this.normalFilename = filename;
    }

    @Option(desc="Tumor BAM file", name="tumor", helpValue="fname")
    public void setTumorFilename(String filename) {
    	this.tumorFilename = filename;
    }
    
    @Option(desc="Regions BED file", name="bed", helpValue="fname")
    public void setBEDFilename(String filename) {
    	this.bedFilename = filename;
    }
    @Option(desc="Region to find copy-number (chr:start-end, using BAM files)", name="region")
    public void setRegion(String region) {
    	this.region = region;
    }
    
    @Option(desc="Normal total read count (optional)", name="norm-total", defaultValue="-1")
    public void setNormalTotal(int count) {
    	this.normalTotal = count;
    }

    @Option(desc="Tumor total read count (optional)", name="tumor-total", defaultValue="-1")
    public void setTumorTotal(int count) {
        	this.tumorTotal = count;
    }

    @Option(desc="Pileup file", name="pileup", helpValue="fname")
    public void setPileupFilename(String pileupFilename) {
        this.pileupFilename = pileupFilename;
    }

	public PileupCopyNumber() {
	}

	public static final double log2(double d) {
		return Math.log(d) / Math.log(2);
	}
	
	@Exec
	public void exec() throws Exception {
		if (pileupFilename == null && (tumorFilename == null || normalFilename == null || (region == null && bedFilename == null))) {
			throw new CommandArgumentException("You must specify either an mpileup file (or stdin), or a normal BAM file, a tumor BAM file, and a region/BED file!");
		}
		
		TabWriter writer = new TabWriter(out);
        writer.write_line("## program: " + CGSeq.getVersion());
        writer.write_line("## cmd: " + CGSeq.getArgs());
		if (pileupFilename != null) {
			writer.write_line("## pileup-input: " + pileupFilename);
		} else {
			writer.write_line("## bed-regions: " + bedFilename);
			writer.write_line("## normal: " + normalFilename);
			writer.write_line("## tumor: " + tumorFilename);
		}
		if (tumorTotal > 0 && normalTotal > 0) {
	        writer.write_line("## tumor-total: " + tumorTotal);
	        writer.write_line("## normal-total: " + normalTotal);
		}
		
		germlineCaller = new SkellamVariantCaller(true, minBaseQual, minMAFDepth);

		if (pileupFilename != null) {
			writer.write("chrom", "start", "end", "ratio (log2)", "copy-number", "norm-median", "tumor-median", "maf_ave", "maf_stdev", "maf_count");
			writer.eol();
			PileupReader reader = new PileupReader(pileupFilename);
			CopyNumberRecord record = calcCopyNumber(reader, normalTotal, tumorTotal);
			if (record != null) {
				writer.write(record.chrom);
				writer.write(record.start);
				writer.write(record.end);
				writer.write(record.ratio);
				writer.write(record.copyNumber);
				writer.write(record.normMedian);
				writer.write(record.tumorMedian);
				writer.write(record.mafAve);
				writer.write(record.mafStdev);
				writer.write(record.mafCount);
				writer.eol();
			}
			reader.close();
		} else if (region == null){
			writer.write("name", "chrom", "start", "end", "ratio (log2)", "copy-number", "norm-median", "tumor-median", "maf_ave", "maf_stdev", "maf_count");
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
				
				ListBuilder<String> lb = new ListBuilder<String>("samtools", "mpileup", "-r", chrom + ":" + ""+(start+1) + "-" + ""+end);
				if (minMapQ > -1) {
					lb.add("-q", ""+minMapQ);
				}
				if (minBaseQual > -1) {
					lb.add("-Q", ""+minBaseQual);
				}
				if (properPairs) {
					lb.add("--rf", "2");
				}
				
				lb.add(normalFilename);
				lb.add(tumorFilename);
				
				ProcessBuilder pb = new ProcessBuilder(lb.list());
				if (verbose) {
					System.err.println(StringUtils.join(" ", pb.command()));
				}
				Process proc = pb.start();
				InputStream bis = new BufferedInputStream(proc.getInputStream());
				PileupReader reader = new PileupReader(bis);

				CopyNumberRecord record = calcCopyNumber(reader, normalTotal, tumorTotal);
				if (record != null) {
					writer.write(name);
					writer.write(chrom);
					writer.write(start);
					writer.write(end);
					writer.write(record.ratio);
					writer.write(record.copyNumber);
					writer.write(record.normMedian);
					writer.write(record.tumorMedian);
					writer.write(record.mafAve);
					writer.write(record.mafStdev);
					writer.write(record.mafCount);
					writer.eol();
				}
				reader.close();
				int retcode = proc.waitFor();
				if (retcode != 0) {
					BufferedReader eis = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
					String line1;
					while ((line1=eis.readLine())!=null) {
						System.err.println("WARNING: " + line1);
					}
					eis.close();
				}
				proc.getErrorStream().close();
			    proc.getInputStream().close();
			    proc.getOutputStream().close();
			    proc.destroy();

			}
			strReader.close();
		} else {
			writer.write("chrom", "start", "end", "ratio (log2)", "copy-number", "norm-median", "tumor-median", "maf_ave", "maf_stdev", "maf_count");
			writer.eol();
			
			ListBuilder<String> lb = new ListBuilder<String>("samtools", "mpileup", "-r", region);
			if (minMapQ > -1) {
				lb.add("-q", ""+minMapQ);
			}
			if (minBaseQual > -1) {
				lb.add("-Q", ""+minBaseQual);
			}
			if (properPairs) {
				lb.add("--rf", "2");
			}

			lb.add(normalFilename);
			lb.add(tumorFilename);
			
			ProcessBuilder pb = new ProcessBuilder(lb.list());
			if (verbose) {
				System.err.println(StringUtils.join(" ", pb.command()));
			}
			
			Process proc = pb.start();
			InputStream bis = new BufferedInputStream(proc.getInputStream());
			PileupReader reader = new PileupReader(bis);

			GenomeSpan gen = GenomeSpan.parse(region, true);
			
			CopyNumberRecord record = calcCopyNumber(reader, normalTotal, tumorTotal);
			if (record != null) {
				writer.write(gen.ref);
				writer.write(gen.start);
				writer.write(gen.end);
				writer.write(tumorTotal);
				writer.write(normalTotal);
				writer.write(record.ratio);
				writer.write(record.copyNumber);
				writer.write(record.normMedian);
				writer.write(record.tumorMedian);
				writer.write(record.mafAve);
				writer.write(record.mafStdev);
				writer.write(record.mafCount);
				writer.eol();
			}
			reader.close();
			int retcode = proc.waitFor();
			if (retcode != 0) {
				BufferedReader eis = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
				String line1;
				while ((line1=eis.readLine())!=null) {
					System.err.println("WARNING: " + line1);
				}
				eis.close();
			}
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
		
		List<Double> maf = new ArrayList<Double>();

		for (PileupRecord pileup: reader) {
			if (chrom == null) {
				chrom = pileup.ref;
				start = pileup.pos;
			}
			
			end = pileup.pos;
			
			normalCounts.add(pileup.getSampleCount(0));
			tumorCounts.add(pileup.getSampleCount(1));
			
			PileupSampleRecord normal = pileup.getSampleRecords(0);
			if (normal.calls == null) {
				continue;
			}

			VariantResults var = germlineCaller.calcVariant(normal.calls, ""); // we don't actually care about REF here.
			
			if (var == null || var.minorCall == null) {
				continue;
			}

			// if there is a minor call, then this is a het.
			// Calculate MAF (not necessarily the B-allele frequency, will always be 0.0-0.5)
			PileupSampleRecord tumor = pileup.getSampleRecords(1);
			
			MapCount<String> counter = new MapCount<String>();
			
			for (PileupBaseCall call: tumor.calls) {
				if (call.qual > minBaseQual || call.op != PileupBaseCallOp.Match) {
					counter.incr(call.toString());
				}
			}
			
			List<Pair<String, Integer>> sorted = counter.getSortedCounts();

			if (sorted.size() < 2) {
				continue;
			}
			
			int major = sorted.get(0).two;
			int minor = sorted.get(1).two;
			
			maf.add(((double) minor / (major + minor)));
		}

		int[] norm = listToArray(normalCounts);
		int[] tumor = listToArray(tumorCounts);
		
		if (norm.length >0 && tumor.length > 0) {
			double normMedian = StatUtils.median(norm);
			double tumorMedian = StatUtils.median(tumor);
			
			double medianRatio = calcCopyRatio(norm, tumor, normalTotal, tumorTotal);
			double copyNumber = calcCopyNumber(medianRatio);

			double mafMean = calcMean(maf);
			double mafStdev = calcStdev(maf, mafMean);
			
			return new CopyNumberRecord(chrom, start, end, medianRatio, copyNumber, normMedian, tumorMedian, mafMean, mafStdev, maf.size());
		}

		return null;
	}
	
	private double calcMean(List<Double> vals) {
		double acc = 0.0;
		double count = 0;
		
		for (Double val: vals) {
			acc += val;
			count++;
		}
		
		return acc / count;
	}

	private double calcStdev(List<Double> vals, double mean) {
		double acc = 0.0;
		double count = 0;
		
		for (Double val: vals) {
			acc += Math.pow((val - mean), 2);
			count++;
		}
		
		double tmp = acc / (count - 1);
		
		return Math.sqrt(tmp);
	}

	
	private int[] listToArray(List<Integer> l) {
		int[] out = new int[l.size()];
		for (int i=0; i<l.size();i++) {
			out[i] = l.get(i);
		}
		return out;
	}
	
	public static double calcCopyNumber(double ratio) {
		return Math.pow(2, ratio + 1);
	}
	public static double calcCopyRatio(int[] normCounts, int[] tumorCounts) {
		return calcCopyRatio(normCounts, tumorCounts, -1, -1);
	}
	
	public static double calcCopyRatio(int[] normCounts, int[] tumorCounts, int normNormalization, int tumorNormalization) {
		double[] normalizedNorm;
		double[] normalizedTumor;

		if (normNormalization > 0 && tumorNormalization > 0) {
			normalizedNorm = normalizeLog2(normCounts, normNormalization);
			normalizedTumor = normalizeLog2(tumorCounts, tumorNormalization);
		} else {
			normalizedNorm = StatUtils.log2(normCounts);
			normalizedTumor = StatUtils.log2(tumorCounts);
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
