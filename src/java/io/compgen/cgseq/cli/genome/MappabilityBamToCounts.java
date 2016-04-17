package io.compgen.cgseq.cli.genome;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.StringLineReader;

import java.util.HashMap;
import java.util.Map;


@Command(name="map-counts", desc="Given a mappability BAM file, determine alignment counts for each region", category="mappability")
public class MappabilityBamToCounts extends AbstractOutputCommand {
	private String filename = "-";
	private String fastaIndex = null;
	
    @UnnamedArg(name = "FILE")
    public void setFilename(String filename) {
    	this.filename = filename;
    }

    @Option(name="fai", desc="Reference FAI file (for chrom sizes)")
    public void setFastaIndex(String fastaIndex) {
    	this.fastaIndex = fastaIndex;
    }

	@Exec
	public void exec() throws Exception {
		Map<String, Integer> refSizes = new HashMap<String, Integer>();
		for (String line: new StringLineReader(fastaIndex)) {
			String[] spl = line.split("\t");
			refSizes.put(spl[0], Integer.parseInt(spl[1]));
		}

		
		
	}
}
