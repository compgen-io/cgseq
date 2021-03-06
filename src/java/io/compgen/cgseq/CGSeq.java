package io.compgen.cgseq;

import io.compgen.cgseq.cli.copynumber.BreakpointFinder;
import io.compgen.cgseq.cli.copynumber.BreakpointMerge;
import io.compgen.cgseq.cli.copynumber.PileupCopyNumber;
import io.compgen.cgseq.cli.genome.Mappability;
import io.compgen.cgseq.cli.genome.MappabilityBamToCounts;
import io.compgen.cgseq.cli.sv.Discord;
import io.compgen.cgseq.cli.varcall.Germline;
import io.compgen.cgseq.cli.varcall.Somatic;
import io.compgen.cgseq.simulation.SimuCall;
import io.compgen.cmdline.Help;
import io.compgen.cmdline.License;
import io.compgen.cmdline.MainBuilder;
import io.compgen.common.StringUtils;

import java.io.IOException;

public class CGSeq {
	public static String getVersion() {
		try {
			return MainBuilder.readFile("io/compgen/cgseq/VERSION");
		} catch (IOException e1) {
			return "unknown";
		}
	}

	private static String args;
	
	public static String getArgs() {
	    return args;
	}

	public static void main(String[] args) throws Exception {
	    CGSeq.args = StringUtils.join(" ", args);
		new MainBuilder()
		.setProgName("cgseq")
		.setHelpHeader("cgseq - Computational Genomics Sequencing Tools\n---------------------------------------")
		.setDefaultUsage("Usage: cgseq cmd [options]")
		.setHelpFooter("http://compgen.io/cgseq\n"+getVersion())
		.addCommand(Help.class)
		.addCommand(License.class)
		.addCommand(PileupCopyNumber.class)
		.addCommand(Discord.class)
		.addCommand(Germline.class)
		.addCommand(Somatic.class)
		.addCommand(SimuCall.class)
		.addCommand(Mappability.class)
		.addCommand(MappabilityBamToCounts.class)
		.addCommand(BreakpointFinder.class)
		.addCommand(BreakpointMerge.class)
		.findAndRun(args);
	}
		
}
