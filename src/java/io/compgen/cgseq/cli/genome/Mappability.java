package io.compgen.cgseq.cli.genome;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import io.compgen.cgseq.CGSUtilsException;
import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.progress.FileChannelStats;
import io.compgen.common.progress.ProgressMessage;
import io.compgen.common.progress.ProgressUtils;

import java.io.File;
import java.io.FileInputStream;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


@Command(name="mappability", desc="Calculate the mappability for a genome (from synthetic BAM file)", category="genome")
public class Mappability extends AbstractOutputCommand {
	
	private class RegionMappabilityScore {
		public final int start;
		public final int end;
		public final int count;
		
		private RegionMappabilityScore(int start, int end, int count) {
			this.start = start;
			this.end = end;
			this.count = count;
		}
		public String toString() {
			return ""+start+"-"+end+"=>"+count;
		}
	}
	
	private class RegionMappabilityScores {
		private List<RegionMappabilityScore> scores = new ArrayList<RegionMappabilityScore>();
		private String currentChrom = null;
		private int currentPos = -1;
		private int currentWriteStart = -1;
		private double currentWriteScore = -1;
		
		private final SamReader reader;
		private double epsilon = 0.00001;
		
		private RegionMappabilityScores(SamReader reader) {
			this.reader = reader;
		}
		
		private void addRegionCount(String chrom, int start, int end, int count) {
//			System.err.println("region-count: "+chrom+":"+start+"-"+end+" => "+count);
			
			if (!chrom.equals(currentChrom)) {
				if (currentChrom != null) {
					clear();
				}
				currentChrom = chrom;
				write(chrom, 0, start, 0);
				currentPos = start;
			}
			this.scores.add(new RegionMappabilityScore(start, end, count));
			calc();
		}

		private void calc() {
			int count=0;
			double acc=0.0;
			
//			System.err.println("calc - currentPos="+currentPos);
			
			for (RegionMappabilityScore score: scores) {
				if (score.start <= currentPos && currentPos < score.end) {
//					System.err.println("    score match: " + score);
					count++;
					acc += score.count;
				}
			}
			
//			System.err.println("  count="+count+", acc="+acc);
			
			if (count > 0) {
				writeScore(currentPos, 1/(acc / count));
			} else {
				writeScore(currentPos, 0);
			}
			
			currentPos++;
			
			while (scores.size() > 0 && scores.get(0).end < currentPos) {
				scores.remove(0);
			}
		}
		
		private void clear() {
			System.err.println("\nclearing - currentChrom:"+currentChrom);
			while (scores.size() > 0) {
				calc();
			}
			
//			if (currentWriteScore > -1) {
//				write(currentChrom, currentWriteStart, currentPos, currentWriteScore);
//				currentWriteStart = currentPos;
//			}

			if (currentPos < reader.getFileHeader().getSequence(currentChrom).getSequenceLength()) {
				write(currentChrom, currentWriteStart, reader.getFileHeader().getSequence(currentChrom).getSequenceLength(), 0);
			}
			
			currentPos = -1;
			currentWriteStart = -1;
			currentWriteScore = -1;

		}
		
		private void writeScore(int pos, double val) {
//			System.err.println("writeScore("+pos+","+val+")");
			double delta = Math.abs(val - currentWriteScore);
			if (delta > epsilon ) {
				if (currentWriteScore > -1) {
					write(currentChrom, currentWriteStart, pos, currentWriteScore);
				}

				currentWriteStart = pos;
				currentWriteScore = val;
			}
		}

		private void write(String chrom, int start, int end, double val) {			
			System.out.println(chrom +"\t"+start +"\t"+end +"\t"+String.format("%.5f", val));
		}
		
	}

	
	private boolean lenient = false;
	private boolean silent = false;
	
	private String filename = "-";
	
	private int maxMismatches = 2;
	private boolean useMAPQ = false;
	private boolean useAS = false;


    @Option(desc="Maximum number of mismatches allowed", name="mismatch", defaultValue="2")
    public void setMaxMismatches(int maxMismatches) {
    	this.maxMismatches = maxMismatches;
    }

    @Option(desc="Use MAPQ values to determine overlap", name="mapq")
    public void setUseMAPQ(boolean useMAPQ) {
    	this.useMAPQ = useMAPQ;
    }

    @Option(desc="Use AS values to determine overlap", name="as")
    public void setUseAS(boolean useAS) {
    	this.useAS = useAS;
    }

    @UnnamedArg(name = "FILE")
    public void setFilename(String filename) throws CommandArgumentException {
    	this.filename = filename;
    }

    @Option(desc = "Use lenient validation strategy", name="lenient")
    public void setLenient(boolean lenient) {
        this.lenient = lenient;
    }

    @Option(desc = "Use silent validation strategy", name="silent")
    public void setSilent(boolean silent) {
        this.silent = silent;
    }    

	@Exec
	public void exec() throws Exception {
		/* in the BAM file, look for all read pairs that map to different regions
		 * region is defined as a chrom:start-end.
		 * 
		 * different regions could be on the same chromosome or different ones.
		 * intragenic SVs need to be at least X bases apart.
		 * 
		 */
	
		// maybe this should be refactored into an iterator of some kind...
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
        if (lenient) {
            readerFactory.validationStringency(ValidationStringency.LENIENT);
        } else if (silent) {
            readerFactory.validationStringency(ValidationStringency.SILENT);
        }

        SamReader reader = null;
        String name;
        FileChannel channel = null;
        if (filename.equals("-")) {
        	throw new CGSUtilsException("BAM input must be from a file, not stdin.");
        } else {
            File f = new File(filename);
            FileInputStream fis = new FileInputStream(f);
            channel = fis.getChannel();
            reader = readerFactory.open(SamInputResource.of(fis));
            name = f.getName();
        }


        Iterator<SAMRecord> it = ProgressUtils.getIterator(name, reader.iterator(), (channel == null)? null : new FileChannelStats(channel), new ProgressMessage<SAMRecord>() {
            @Override
            
            public String msg(SAMRecord current) {
            	return current.getReadName();
            }});

        RegionMappabilityScores scores = new RegionMappabilityScores(reader);
        
        String currentName = null;
        int count = 0;
        int bufVal = -1;
        
        while (it.hasNext()) {
            SAMRecord read = it.next();
            if (read.getReadUnmappedFlag()) {
            	// shouldn't happen
            	continue;
            }

            
            if (!read.getReadName().equals(currentName)) {
            	if (currentName != null) {
                    String[] spl = currentName.split(":");
            		String[] se = spl[1].split("-");
            		
            		String chrom = spl[0];
            		int start = Integer.parseInt(se[0]);
            		int end = Integer.parseInt(se[1]);

            		scores.addRegionCount(chrom, start, end, count);
            	}
            	
            	currentName = read.getReadName();
            	count = 0;
            }
            
            if (useMAPQ) {
            	if (bufVal == -1 || read.getMappingQuality() > bufVal) {
            		bufVal = read.getMappingQuality();
            		count = 1;
            	} else if (read.getMappingQuality() == bufVal) {
            		count ++;
            	}
            } else if (useAS) {
            	int asVal = read.getIntegerAttribute("AS");
            	if (bufVal == -1 || asVal > bufVal) {
            		bufVal = asVal;
            		count = 1;
            	} else if (asVal == bufVal) {
            		count ++;
            	}
            } else {
            	int mismatches = read.getIntegerAttribute("NM");

            	for (CigarElement el:read.getCigar().getCigarElements()) {
            		if (el.getOperator() != CigarOperator.M) {
            			mismatches += el.getLength();
            		}
            	}
            	
//            	System.err.println(read.getReadName() + "; " + read.getReferenceName()+":"+read.getAlignmentStart()+" => nm="+ read.getIntegerAttribute("NM") + ", cigar="+read.getCigarString()+", mismatches="+mismatches);
            	
            	
            	if (mismatches <= maxMismatches) {
            		count++;
            	}
            }
            
        }

    	if (currentName != null) {
    		String[] spl = currentName.split(":");
    		String[] se = spl[1].split("-");
    	
    		scores.addRegionCount(spl[0], Integer.parseInt(se[0]), Integer.parseInt(se[1]), count);
    	}
    	
    	scores.clear();
    	

        reader.close();

	}
}
