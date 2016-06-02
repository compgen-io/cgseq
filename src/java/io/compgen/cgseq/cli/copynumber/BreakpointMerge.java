package io.compgen.cgseq.cli.copynumber;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.ComparablePair;
import io.compgen.common.StringLineReader;
import io.compgen.common.StringUtils;
import io.compgen.common.TabWriter;
import io.compgen.ngsutils.support.stats.StatUtils;
import io.compgen.ngsutils.support.stats.StatUtils.MeanStdDev;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

@Command(name="bp-calc", desc="Given bp-dist file(s), calculate potential breakpoints", category="copy-number")
public class BreakpointMerge extends AbstractOutputCommand {
	private String[] filenames = null;
	private String gapFile = null;
	private double sigmas = 2.0;
	private double usePct = .95;

	
	@Option(desc="Threshold for bp-dist (in sigmas)", charName="s", defaultValue="2.0")
    public void setRegion(double sigma) {
    	this.sigmas = sigma;
    }

	@Option(name="pct", desc="User the bottom percentage of sites for mean/threshold calc (0-1.0)", defaultValue="0.95")
    public void setUsePct(double usePct) {
    	this.usePct = usePct;
    }

	@Option(desc="Genome gaps to include (BED)", name="gaps")
	public void setGapFile(String gapFile) {
		this.gapFile = gapFile;
	}

    @UnnamedArg(name = "FILE1 FILE2 ...")
    public void setFilename(String[] filenames) throws CommandArgumentException {
        this.filenames = filenames;
    }


	@Exec
	public void exec() throws CommandArgumentException, IOException {
		if (filenames == null) {
			throw new CommandArgumentException("You must specify at least one bp-dist file.");
		}
		
		Map<ComparablePair<String, Integer>, Double> bpDists = new TreeMap<ComparablePair<String, Integer>, Double>();
		Map<String, Integer> seqLengths = new HashMap<String, Integer>();
		
		for (String fname: filenames) {
			System.err.println("Reading: " + fname);
			StringLineReader reader = new StringLineReader(fname);
			for (String line: reader) {
				if (line == null) {
					continue;
				}
				
				if (line.charAt(0) == '#') {
					if (line.startsWith("## ref ")) {
						String[] spl = line.split(" ");
						if (!seqLengths.containsKey(spl[2])) {
							seqLengths.put(spl[2], Integer.parseInt(spl[3]));
						}
					}
				}
				
				String[] cols = StringUtils.strip(line).split("\t");
				if (cols.length>4 && !cols[3].equals("")) {
					String chrom = cols[0];
					int pos = Integer.parseInt(cols[4]);
					double dist = Double.parseDouble(cols[3]);
					
					ComparablePair<String, Integer> chromPos = new ComparablePair<String, Integer>(chrom, pos);
					if (bpDists.containsKey(chromPos)) {
						if (bpDists.get(chromPos) < dist) {
							bpDists.put(chromPos,  dist);
						}
					} else {
						bpDists.put(chromPos,  dist);
					}
				}
			}
		}
				
		
		double[] allDistances = new double[bpDists.size()];
		
		{
			int i=0;
			for (Double dist: bpDists.values()) {
				allDistances[i++]=dist;
			}
		}

		Arrays.sort(allDistances);
		double[] distances = new double[(int) (allDistances.length * usePct)];
		for (int i=0; i<distances.length; i++) {
			distances[i] = allDistances[i];
		}


		MeanStdDev msdAll = StatUtils.calcMeanStdDev(allDistances);
		
		System.err.println("Mean    (all) : " + msdAll.mean);
		System.err.println("Std-dev (all) : " + msdAll.stddev);
		
		MeanStdDev msd = StatUtils.calcMeanStdDev(distances);
		double threshold = msd.mean + (sigmas * msd.stddev);
		
		System.err.println("Mean     : " + msd.mean);
		System.err.println("Std-dev  : " + msd.stddev);
		System.err.println("Threshold: " + threshold);

		TabWriter writer = new TabWriter(out);
		
		int curPos = 0;
		String curChrom = null;
		int count = 0;
		
		Set<String> usedChroms = new HashSet<String>();
		
		for (Entry<ComparablePair<String, Integer>, Double> e: bpDists.entrySet()) {
			if (e.getValue() > threshold) {
				String chrom = e.getKey().one;
				Integer newPos = e.getKey().two;
				
				if (curChrom == null || !curChrom.equals(chrom)) {
					if (curChrom != null) {
						writer.write(curChrom);
						writer.write(curPos);
						writer.write(seqLengths.get(curChrom));
						writer.write("region_"+(++count));
						writer.write("0");
						writer.eol();
					}
					curChrom = chrom;
					usedChroms.add(chrom);
					curPos = 0;
				}

				writer.write(curChrom);
				writer.write(curPos);
				writer.write(newPos);
				writer.write("region_"+(++count));
				writer.write(e.getValue());
				writer.eol();
				curPos = newPos;
			}
		}
		writer.write(curChrom);
		writer.write(curPos);
		writer.write(seqLengths.get(curChrom));
		writer.write("region_"+(++count));
		writer.write("0");
		writer.eol();
		
		for (String chrom: seqLengths.keySet()) {
			if (!usedChroms.contains(chrom)) {
				writer.write(chrom);
				writer.write(0);
				writer.write(seqLengths.get(chrom));
				writer.write("region_"+(++count));
				writer.write("0");
				writer.eol();
			}
		}
		
		writer.close();
		System.err.println("Total regions: " + count);
	}
}
