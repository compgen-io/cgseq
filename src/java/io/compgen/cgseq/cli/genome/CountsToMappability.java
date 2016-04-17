package io.compgen.cgseq.cli.genome;

import htsjdk.samtools.CigarElement;
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
import java.util.Iterator;


@Command(name="map-counts", desc="Given a mappability BAM file, determine alignment counts for each region", category="mappability")
public class CountsToMappability extends AbstractOutputCommand {
	private boolean lenient = false;
	private boolean silent = false;
	
	private String filename = "-";
	
	private int maxMismatches = 2;
	private boolean useMAPQ = false; // default is to use mismatches...
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
            		System.out.println(spl[0]+"\t" + se[0]+ "\t" + se[1] + "\t" + count);
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
            		switch (el.getOperator()) {
            		case D:
            		case I:
            		case S:
            		case N:
            			mismatches += el.getLength();
					default:
						break;
            		}
            	}
            	
            	if (mismatches <= maxMismatches) {
            		count++;
            	}
            }
            
        }

    	if (currentName != null) {
    		String[] spl = currentName.split(":");
    		String[] se = spl[1].split("-");
    		System.out.println(spl[0]+"\t" + se[0]+ "\t" + se[1] + "\t" + count);
    	}
    	

        reader.close();

	}
}
