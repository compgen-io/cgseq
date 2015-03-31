package io.compgen.cgseq.cli.varcall;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import io.compgen.cgseq.CGSeq;
import io.compgen.cgseq.variant.PoissonVariant;
import io.compgen.cgseq.variant.VariantResults;
import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.IterUtils;
import io.compgen.common.StringUtils;
import io.compgen.common.TabWriter;
import io.compgen.ngsutils.pileup.BAMPileup;
import io.compgen.ngsutils.pileup.PileupRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


@Command(name="germline", desc="Call variants for a germline sample (diploid)", category="variants")
public class Germline extends AbstractOutputCommand {
	private String filename = "-";
	private String refFilename;
	
	private int minBaseQual = 30;
	private int minMappingQual = 0;

	private int filterFlags = -1;
    private int requiredFlags = -1;

    private boolean disableBAQ = true;
    private boolean extendedBAQ = false;
    private boolean onlyVariants = false;
    
    private double expectedAlleleFrequency = 0.5;
    private boolean backgroundCorrect = true;
    private int minCoverage = 10;
    
    private double minMinorStrandFreq = 0.05;
    
    @Option(desc="Filter flags", name="filter-flags")
    public void setFilterFlags(int filterFlags) {
    	this.filterFlags = filterFlags;
    }
    
    @Option(desc="Required flags", name="required-flags")
    public void setRequiredFlags(int requiredFlags) {
    	this.requiredFlags = requiredFlags;
    }
    
    @Option(desc="Only show variants", name="only-variants")
    public void setOnlyVariants(boolean val) {
    	this.onlyVariants = val;
    }
    
    @Option(desc="No background correction", name="nobg")
    public void setNoBackgroundCorrect(boolean val) {
    	this.backgroundCorrect = val;
    }

    @Option(desc="Minimum coverage", name="mincoverage", defaultValue="10")
    public void setMinCoverage(int minCoverage) {
    	this.minCoverage = minCoverage;
    }

    @Option(desc="Estimated allele frequency", name="allelefreq", defaultValue="0.5")
    public void setEstAlleleFreq(double estAlleleFreq) {
    	this.expectedAlleleFrequency = estAlleleFreq;
    }

    @Option(desc="Minimum base-quality score", name="basequal", defaultValue="30")
    public void setMinBaseQual(int minBaseQual) {
    	this.minBaseQual = minBaseQual;
    }

    @Option(desc="Minimum alignment mapping score", name="mapqual", defaultValue="0")
    public void setMinMapQual(int minMappingQual) {
    	this.minMappingQual = minMappingQual;
    }

    @Option(desc = "Skip BAQ recalculation", name="disable-baq")
    public void setDisableBAQ(boolean val) {
    	this.disableBAQ = val;
    }

    @Option(desc = "Extended BAQ calculation", name="extended-baq")
    public void setExtendedBAQ(boolean val) {
    	this.extendedBAQ = val;
    }
    
    @UnnamedArg(name = "ref_fasta germline_bam", required=true)
    public void setFilename(List<String> filenames) throws CommandArgumentException {
    	if (filenames.size()!=2) {
            throw new CommandArgumentException("You must specify a reference genome (FASTA) and germline BAM file!");
    	}
        this.refFilename = filenames.get(0);
    	this.filename = filenames.get(1);
    }

	public Germline() {
	}

	@Exec
	public void exec() throws Exception {
//				##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
//				##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
//				##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
//				##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
//				##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
//				##INFO=<ID=IS,Number=2,Type=Float,Description="Maximum number of reads supporting an indel and fraction of indel reads">
//				##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
//				##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
//				##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
//				##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
//				##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
//				##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
//				##INFO=<ID=PV4,Number=4,Type=Float,Description=	"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
//				##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
//				##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
//				##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
//				##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
//				##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
//				##INFO=<ID=QBD,Number=1,Type=Float,Description="Quality by Depth: QUAL/#reads">
//				##INFO=<ID=RPB,Number=1,Type=Float,Description="Read Position Bias">
//				##INFO=<ID=MDV,Number=1,Type=Integer,Description="Maximum number of high-quality nonRef reads in samples">
//				##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias (v2) for filtering splice-site artefacts in RNA-seq data. Note: this version may be broken
//				##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
//				##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
//				##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
//				##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
//				##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality non-reference bases">
//				##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
//				##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">


		BAMPileup pileup = new BAMPileup(filename);
		pileup.setDisableBAQ(disableBAQ);
		pileup.setExtendedBAQ(extendedBAQ);
		pileup.setFlagFilter(filterFlags);
		pileup.setFlagRequired(requiredFlags);
		pileup.setMinBaseQual(minBaseQual);
		pileup.setMinMappingQual(minMappingQual);
		pileup.setRefFilename(refFilename);

		TabWriter writer = new TabWriter(out);
		writer.write_line("##fileformat=VCFv4.1");
		writer.write_line("##cgseqVersion="+CGSeq.getVersion());
		writer.write_line("##reference=file://"+new File(refFilename).getAbsolutePath());
		writer.write_line("##pileupCommand="+StringUtils.join(" ", pileup.getCommand()));

		SamReader bam = SamReaderFactory.makeDefault().open(new File(filename));
		SAMFileHeader header = bam.getFileHeader();
		for (SAMSequenceRecord rec: header.getSequenceDictionary().getSequences()) {
			writer.write_line("##contig=<ID="+rec.getSequenceName()+",length="+rec.getSequenceLength()+">");
		}
		bam.close();

		writer.write_line("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
		writer.write_line("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">");

		writer.write("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", filename);
		writer.eol();

		PoissonVariant poissonCaller = new PoissonVariant(backgroundCorrect, minBaseQual, expectedAlleleFrequency);
		
		for (PileupRecord record: IterUtils.wrap(pileup.pileup())) {
			if (record.getSampleRecords(0).coverage < minCoverage) {
				continue;
			}
			
			// Germline should be the first record
			VariantResults varResult = poissonCaller.calcVariant(record.getSampleRecords(0).calls);
			
			if (varResult == null) {
				continue;
			}

			// Assume the following possible genotypes: REF:REF, REF:ALT, ALT:ALT.

			String altCall = null;
			int qual=-1;
			int DP=-1;
			String DP4=""; 

			// Is this a REF:ALT call?
			if (varResult.altPvalue > varResult.majorPvalue) {
				if (varResult.getMinorStrandFreq() >= minMinorStrandFreq) {
					// Yes - AB
					qual = pvalueToPhred(varResult.altError);
					DP = varResult.majorCount + varResult.minorCount;
					if (record.refBase.toUpperCase().equals(varResult.majorCall.toUpperCase())) {
						// REF:ALT
						altCall = varResult.minorCall;
						DP4 = varResult.majorPlusStrandCount+","+varResult.majorMinusStrandCount+","+varResult.minorPlusStrandCount+","+varResult.minorMinusStrandCount;
					} else {
						// ALT:REF
						altCall = varResult.majorCall;
						DP4 = varResult.minorPlusStrandCount+","+varResult.minorMinusStrandCount+","+varResult.majorPlusStrandCount+","+varResult.majorMinusStrandCount;
					}
				}
			}
			
			if (altCall == null) {
				DP = varResult.majorCount;
				qual = pvalueToPhred(varResult.majorError);
				if (!record.refBase.toUpperCase().equals(varResult.majorCall.toUpperCase())) {
					// TODO: check for multi-allele? 
					if (varResult.getMajorStrandFreq() >= minMinorStrandFreq) {
						altCall = varResult.majorCall;
						DP4 = "0,0,"+varResult.majorPlusStrandCount+","+varResult.majorMinusStrandCount; // ALT:ALT
					}
				}
				
				if (altCall == null) {
					DP4 = varResult.majorPlusStrandCount+","+varResult.majorMinusStrandCount+",0,0"; // REF:REF
				}
			}

			// otherwise we could be REF:REF or ALT:ALT.

			if (onlyVariants && altCall == null) {
				continue;
			}
			
			
			if (altCall != null && altCall.equals(varResult.majorCall)) {
				
			}

			
			
			
			
			
			writer.write(record.ref);
			writer.write(record.pos+1);
			writer.write("."); // dbsnp id
			writer.write(record.refBase);
			
			if (altCall != null) {
				writer.write(altCall);
			} else {
				writer.write(".");
			}

			// qual is prob we are wrong (for either way...)
			writer.write(qual);
			
			// info
			List<String> info = new ArrayList<String>();
			info.add("DP="+DP);
			info.add("DP4="+DP4);
			info.add("CALLS="+varResult.majorCall+","+varResult.minorCall);
			info.add("COUNTS="+varResult.majorCount+","+varResult.minorCount);
			info.add("MAJSTRAND_FREQ="+varResult.getMajorStrandFreq());
			info.add("MINSTRAND_FREQ="+varResult.getMinorStrandFreq());
			
			writer.write(StringUtils.join(";", info));
			
			
			// format
			writer.write("");
//			writer.write(varResult.majorCall);
//			writer.write(varResult.minorCall);
//			writer.write(varResult.majorCount);
//			writer.write(varResult.minorCount);
//			writer.write(varResult.bgCount);
//			writer.write(varResult.hetPvalue);
//			writer.write(varResult.homPvalue);
//			writer.write(varResult.majorPlusStrandFreq);
//			writer.write(varResult.minorPlusStrandFreq);
//			writer.write(varResult.majorPlusPvalue);
//			writer.write(varResult.minorPlusPvalue);
			writer.eol();
		}
		writer.close();
	}
	
	public int pvalueToPhred(double pval) {
		return (int) (-10 * Math.log10(pval));
	}
}
