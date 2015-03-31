package io.compgen.cgseq;

import io.compgen.cgseq.cli.copynumber.Breakpoints;
import io.compgen.cgseq.cli.copynumber.PileupCopyNumber;
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
			return MainBuilder.readFile("VERSION");
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
		.setHelpHeader("cgseq - Computational Genomics Sequencing Analysis\n---------------------------------------")
		.setDefaultUsage("Usage: cgseq cmd [options]")
		.setHelpFooter("http://compgen.io/cgseq\n"+MainBuilder.readFile("VERSION"))
		.setCategoryOrder(new String[]{"server", "client", "help"})
		.addCommand(Help.class)
		.addCommand(License.class)
		.addCommand(Breakpoints.class)
		.addCommand(PileupCopyNumber.class)
		.addCommand(Discord.class)
		.addCommand(Germline.class)
		.addCommand(Somatic.class)
		.addCommand(SimuCall.class)
		.findAndRun(args);
	}
		
}
	
//	static private Map<String, Class<Exec>> execs = new HashMap<String, Class<Exec>>();
//	static {
//		loadExec(Breakpoints.class);
//		loadExec(PileupCopyNumber.class);
//	}
//
//	@SuppressWarnings("unchecked")
//	private static void loadExec(Class<?> cls) {
//		String name = cls.getName().toLowerCase();
//		Command cmd = (Command) cls.getAnnotation(Command.class);
//		if (cmd != null) {
//			name = cmd.name();
//		}
//		execs.put(name, (Class<Exec>) cls);
//	}
//
//	public static void usage() {
//		usage(null);
//	}
//
//	public static void usage(String msg) {
//		if (msg != null) {
//			System.err.println(msg);
//			System.err.println();
//		}
//		System.err.println("CGSeq - Tools for cancer genome analysis");
//		System.err.println("");
//		System.err.println("Usage: cgsutils cmd {options}");
//		System.err.println("");
//		System.err.println("Available commands:");
//		int minsize = 12;
//		String spacer = "            ";
//		for (String cmd : execs.keySet()) {
//			if (cmd.length() > minsize) {
//				minsize = cmd.length();
//			}
//		}
//		Map<String, List<String>> progs = new HashMap<String, List<String>>();
//
//		for (String cmd : execs.keySet()) {
//			Command c = execs.get(cmd).getAnnotation(Command.class);
//			if (c != null) {
//				if (!progs.containsKey(c.cat())) {
//					progs.put(c.cat(), new ArrayList<String>());
//				}
//
//				if (!c.desc().equals("")) {
//					spacer = "";
//					for (int i = cmd.length(); i < minsize; i++) {
//						spacer += " ";
//					}
//					spacer += " - ";
//					progs.get(c.cat()).add("  " + cmd + spacer + c.desc());
//				} else {
//					progs.get(c.cat()).add("  " + cmd);
//				}
//			} else {
//				if (!progs.containsKey("General")) {
//					progs.put("General", new ArrayList<String>());
//				}
//				progs.get("General").add("  " + cmd);
//
//			}
//		}
//
//		List<String> cats = new ArrayList<String>(progs.keySet());
//		Collections.sort(cats);
//
//		for (String cat : cats) {
//            System.err.println("[" + cat + "]");
//			Collections.sort(progs.get(cat));
//			for (String line : progs.get(cat)) {
//				System.err.println(line);
//			}
//            System.err.println("");
//		}
//
//		spacer = "";
//		for (int i = 12; i < minsize; i++) {
//			spacer += " ";
//		}
//		spacer += " - ";
//		System.err.println("[help]");
//		System.err.println("  help command" + spacer
//				+ "Help message for the given command");
//		
//		System.err.println("");
//		System.err.println(getVersion());
//	}
//	
//	
//	public static void main(String[] args) throws Exception {
//	    
//		if (args.length == 0) {
//			usage();
//		} else if (args[0].equals("help")) {
//			if (args.length == 1) {
//				usage();
//			} else {
//				if (!execs.containsKey(args[1])) {
//					usage("Unknown command: " + args[1]);
//				} else {
//					showHelp(execs.get(args[1]));
//				}
//			}
//		} else if (execs.containsKey(args[0])) {
//			List<String> l = Arrays.asList(args).subList(1, args.length);
//			try {
//				Exec exec = CliFactory.parseArgumentsUsingInstance(execs
//						.get(args[0]).newInstance(), (String[]) l
//						.toArray(new String[l.size()]));
//				exec.exec();
//			} catch (HelpRequestedException e) {
//				System.err.println(e.getMessage());
//				System.err.println("");
//				System.err.println(getVersion());
//			} catch (ArgumentValidationException e) {
//				System.err.println(e.getMessage());
//				showHelp(execs.get(args[0]));
//			} catch (Throwable t) {
//				t.printStackTrace(System.err);
//				System.exit(1);
//			}
//		} else {
//			usage("Unknown command: " + args[0]);
//		}
//	}
//
//	private static void showHelp(Class<Exec> clazz) throws Exception {
//		Command cmd = clazz.getAnnotation(Command.class);
//		if (cmd != null) {
//			if (cmd.desc().equals("")) {
//				System.err.println(cmd.name());
//			} else {
//				System.err.println(cmd.name() + " - " + cmd.desc());
//			}
//			System.err.println("");
//
//			if (!cmd.doc().equals("")) {
//				System.err.println(cmd.doc());
//			}
//		} else {
//			System.err.println(clazz.getName().toLowerCase());
//		}
//		String[] args = { "--help" };
//		try {
//			CliFactory.parseArgumentsUsingInstance(clazz.newInstance(), args);
//		} catch (HelpRequestedException e) {
//			System.err.println(e.getMessage());
//		}
//		System.err.println("");
//		System.err.println(getVersion());
//	}
//}
