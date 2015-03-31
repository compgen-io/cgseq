package io.compgen.cgseq.cli.varcall;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;


@Command(name="somatic", desc="Call variants for a somatic sample (vs germline)", category="variants")
public class Somatic extends AbstractOutputCommand {
	private String filename = "-";
	private String germlineFilename = "-";
	private String refFilename;
	
	private int minBaseQual = 30;
	private int minMappingQual = 0;

	private int filterFlags = -1;
    private int requiredFlags = -1;

    private boolean disableBAQ = true;
    private boolean extendedBAQ = false;

    @Option(desc="Filter flags", name="filter-flags")
    public void setFilterFlags(int filterFlags) {
    	this.filterFlags = filterFlags;
    }
    
    @Option(desc="Required flags", name="required-flags")
    public void setRequiredFlags(int requiredFlags) {
    	this.requiredFlags = requiredFlags;
    }
    
    @Option(desc="Minimum base-quality score", name="basequal", defaultValue="13")
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
    
    @UnnamedArg(name = "ref_fasta germline_bam somatic_bam")
    public void setFilename(String[] filenames) throws CommandArgumentException {
    	if (filenames.length!=3) {
            throw new CommandArgumentException("You must specify a reference genome (FASTA), germline (BAM), and somatic (BAM) files!");
    	}
        this.refFilename = filenames[0];
    	this.germlineFilename = filenames[1];
    	this.filename = filenames[2];
    }

	@Exec
	public void exec() throws Exception {
		
	}
}
