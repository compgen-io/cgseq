package io.compgen.cgseq.cli.sv;

import io.compgen.ngsutils.annotation.GenomeSpan;

import java.util.HashSet;
import java.util.Set;

public class Translocation {
	private GenomeSpan from;
	private GenomeSpan to;
	private Set<String> readnames = new HashSet<String>();
	private int posStrandCount = 0;
	private int negStrandCount = 0;
	private Translocation reciprocal = null;
	
	public Translocation(GenomeSpan from, GenomeSpan to) {
		this.from = from;
		this.to = to;
	}
	
	public void addEvidence(String readName, boolean isNegStrand) {
		readnames.add(readName);
		if (isNegStrand) {
			negStrandCount++;
		} else {
			posStrandCount++;
		}
	}
	
	public int getEvidenceCount() {
		return readnames.size();
	}
	
	public int getNegStrandCount() {
		return negStrandCount;
	}
	
	public int getPosStrandCount() {
		return posStrandCount;
	}
	
	public boolean containsRead(String readName) {
		return readnames.contains(readName);
	}
	
	public boolean overlapsFrom(GenomeSpan span) {
		return from.overlaps(span); 
	}

	public boolean overlapsTo(GenomeSpan span) {
		return to.overlaps(span); 
	}

	public void extendFrom(GenomeSpan span) {
		int start = Math.min(from.start, span.start);
		int end = Math.max(from.end,  span.end);
		from = new GenomeSpan(from.ref, start, end);
	}

	public void extendTo(GenomeSpan span) {
		int start = Math.min(to.start, span.start);
		int end = Math.max(to.end,  span.end);
		to = new GenomeSpan(to.ref, start, end);
	}
	
	public GenomeSpan getFrom() {
		return from;
	}
	
	public GenomeSpan getTo() {
		return to;
	}
	
	public String toString() {
		return from.ref+"\t"+from.start+"\t"+from.end+"\t"+to.ref+"\t"+to.start+"\t"+to.end + "\t"+posStrandCount+"\t"+negStrandCount;
	}

	public void combine(Translocation trans) {
		int fstart = Math.min(trans.from.start, from.start);
		int fend = Math.max(trans.from.end, from.end);
		from = new GenomeSpan(from.ref, fstart, fend);
	
		int tstart = Math.min(trans.to.start, to.start);
		int tend = Math.max(trans.to.end, to.end);
		to = new GenomeSpan(to.ref, tstart, tend);
		
		negStrandCount += trans.negStrandCount;
		posStrandCount += trans.posStrandCount;
		
		readnames.addAll(trans.readnames);
	}
	
	public Translocation getReciprocal() {
		return this.reciprocal;
	}
	public void setReciprocal(Translocation trans) {
		this.reciprocal = trans;
	}
}
