package io.compgen.cgseq.cli.sv;

import htsjdk.samtools.SAMRecord;
import io.compgen.ngsutils.annotation.GenomeSpan;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Translocations {
	public class Translocation {
		private GenomeSpan from;
		private GenomeSpan to;
		private Set<String> readnames = new HashSet<String>();
		
		public Translocation(GenomeSpan from, GenomeSpan to) {
			this.from = from;
			this.to = to;
		}
		
		public void addEvidence(String readName) {
			readnames.add(readName);
		}
		
		public boolean containsRead(String readName) {
			return readnames.contains(readName);
		}
		
		public boolean overlaps(GenomeSpan span) {
			return from.overlaps(span) || to.overlaps(span); 
		}

		public void extend(GenomeSpan span) {
			if (from.overlaps(span)) {
				int start = Math.min(from.start, span.start);
				int end = Math.max(from.end,  span.end);
				from = new GenomeSpan(from.ref, start, end, from.strand);
			} else if (to.overlaps(span)) {
				int start = Math.min(to.start, span.start);
				int end = Math.max(to.end,  span.end);
				to = new GenomeSpan(to.ref, start, end, to.strand);
			}
			
		}
		
		public String toString() {
			return from.ref+":"+from.start+"-"+from.end+" => "+to.ref+":"+to.start+"-"+to.end + " ("+readnames.size()+")";
		}
	}
	
	private List<Translocation> translocations = new ArrayList<Translocation>();
	private final int extend;
	
	public Translocations(int extend) {
		this.extend = extend;
	}

	public void addRead(SAMRecord read) {
    	GenomeSpan from = new GenomeSpan(read.getReferenceName(), read.getAlignmentStart()-1-extend, read.getAlignmentEnd()+extend);

    	// the "to" is only an estimate of the proper span for now... when we parse the mate (later in the file), this will get fixed.
    	GenomeSpan to = new GenomeSpan(read.getMateReferenceName(), read.getMateAlignmentStart()-1-extend, read.getMateAlignmentStart()+extend);
    
    	addTranslocation(read.getReadName(), from, to);
	}

	public void addTranslocation(String readName, GenomeSpan from, GenomeSpan to) {
		// TODO: toss translocations that are out of this window?
		
		
		for (Translocation t: translocations) {
			if (t.containsRead(readName)) {
				t.extend(from);
				t.extend(to);
				t.addEvidence(readName);
				// found by name, return
				return;
			}
		}

		for (Translocation t: translocations) {
			if (t.overlaps(from) || t.overlaps(to)) {
				t.extend(from);
				t.extend(to);
				t.addEvidence(readName);
				return;
			}
		}
		
		Translocation t = new Translocation(from, to);
		t.addEvidence(readName);
		translocations.add(t);
	}
	
	public int size() {
		return translocations.size();
	}
	
	public void dump() {
		for (Translocation t: translocations) {
			System.out.println(t);
		}
	}
}