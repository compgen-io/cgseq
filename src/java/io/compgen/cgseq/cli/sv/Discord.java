package io.compgen.cgseq.cli.sv;

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
import io.compgen.common.Counter;
import io.compgen.common.IterUtils;
import io.compgen.common.IterUtils.MapFunc;
import io.compgen.common.StringUtils;
import io.compgen.common.progress.FileChannelStats;
import io.compgen.common.progress.ProgressMessage;
import io.compgen.common.progress.ProgressUtils;
import io.compgen.ngsutils.bam.support.ReadUtils;

import java.io.File;
import java.io.FileInputStream;
import java.nio.channels.FileChannel;
import java.util.Iterator;


@Command(name="discord", desc="Extract discordant reads from a BAM file", category="sv")
public class Discord extends AbstractOutputCommand {
	private boolean lenient = false;
	private boolean silent = false;
	
	private String filename = "-";
	private String refFilename = null;
	
	private int minMappingQual = 0;

	private int filterFlags = -1;
    private int requiredFlags = -1;

    private int minIntraChromDist = 5000; // for DNA this is probably okay. for RNA, this should be closer to 200_000 - 500_000.
    private int extendBases = 1000;
    private int minReadCount = 2;
    
    @Option(desc="Filter flags", name="filter-flags")
    public void setFilterFlags(int filterFlags) {
    	this.filterFlags = filterFlags;
    }
    
    @Option(desc="Required flags", name="required-flags")
    public void setRequiredFlags(int requiredFlags) {
    	this.requiredFlags = requiredFlags;
    }
    
    @Option(desc="Minimum intrachromasomal distance", name="intradist", defaultValue="5000")
    public void setMinIntraChromDist(int minIntraChromDist) {
    	this.minIntraChromDist = minIntraChromDist;
    }

    @Option(desc="Extend target bases (for finding reciprocal regions)", name="extend-bases", defaultValue="1000")
    public void setExtendBases(int extendBases) {
    	this.extendBases = extendBases;
    }

    @Option(desc="Minimum alignment mapping score", name="mapqual", defaultValue="0")
    public void setMinMapQual(int minMappingQual) {
    	this.minMappingQual = minMappingQual;
    }

    @Option(desc="Minimum number of supporting reads", name="minreads", defaultValue="2")
    public void setMinReadCount(int minReadCount) {
    	this.minReadCount = minReadCount;
    }

    @UnnamedArg(name = "ref_fasta mapped_bam")
    public void setFilename(String[] filenames) throws CommandArgumentException {
    	if (filenames.length!=2) {
            throw new CommandArgumentException("You must specify a reference genome (FASTA) and mapped reads (BAM) files!");
    	}
        this.refFilename = filenames[0];
    	this.filename = filenames[1];
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

        final Translocations translocations = new Translocations(extendBases, minReadCount);
        final Counter total = new Counter();

        Iterator<SAMRecord> it = ProgressUtils.getIterator(name, reader.iterator(), (channel == null)? null : new FileChannelStats(channel), new ProgressMessage<SAMRecord>() {
            @Override
            public String msg(SAMRecord current) {
                return total.getValue() + " (" + translocations.size()+") "+current.getReferenceName()+":"+current.getAlignmentStart();
            }});

        while (it.hasNext()) {
        	total.incr();
            
            SAMRecord read = it.next();
            if (!read.getReadPairedFlag() || read.getReadUnmappedFlag() || read.getMateUnmappedFlag() || read.getDuplicateReadFlag()) {
            	// We require paired end data for this...
            	// For the first round, both pairs need to have mapped.
            	continue;
            }
            
            if (!ReadUtils.isReadUniquelyMapped(read)) {
            	// If this is a repetitive read, skip it... too many false positives. 
            	continue;
            }
            
            if (requiredFlags > 0 && ((read.getFlags() & requiredFlags) != requiredFlags)) {
            	// missing a required flag
            	continue;
            }

            if (filterFlags > 0 && ((read.getFlags() & filterFlags) > 0)) {
            	// has a bad flag
            	continue;
            }

            if (this.minMappingQual > 0 && read.getMappingQuality() < this.minMappingQual) {
            	continue;
            }
            
            if (read.getReferenceIndex() != read.getMateReferenceIndex() || read.getMateAlignmentStart() - read.getAlignmentEnd() > minIntraChromDist) {
            	// discordant read!
            	translocations.addRead(read);
            }
//            if (total.getValue() % 10000 == 0) {
//            	System.err.println();
//            	System.err.println("=================================");
//	            System.err.println(StringUtils.join("\n", IterUtils.map(translocations.peekBuffer(), new MapFunc<Translocation, String>() {
//					@Override
//					public String map(Translocation obj) {
//						return obj.getFrom().toString()+"=>"+obj.getTo().toString();
//					}})));
//            }
        }
        translocations.done();
        reader.close();
        System.err.println("Successfully read: "+total.getValue()+" records.");

        translocations.findReciprocalMatches();

        translocations.dump();		
	}
}
