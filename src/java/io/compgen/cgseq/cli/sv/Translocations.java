package io.compgen.cgseq.cli.sv;

import htsjdk.samtools.SAMRecord;
import io.compgen.common.IterUtils;
import io.compgen.common.progress.IncrementingStats;
import io.compgen.common.progress.ProgressUtils;
import io.compgen.ngsutils.annotation.AbstractAnnotationSource;
import io.compgen.ngsutils.annotation.GenomeAnnotation;
import io.compgen.ngsutils.annotation.GenomeSpan;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Translocations extends AbstractAnnotationSource<Translocation> {
	@Override
	public String[] getAnnotationNames() {
		return null;
	}
	
	private final int extend;
	private final int minEvidence;
	private List<Translocation> buffer = new ArrayList<Translocation>();
	
	public Translocations(int extend, int minEvidence) {
		this.extend = extend;
		this.minEvidence = minEvidence;
	}

	public void addRead(SAMRecord read) {
		// Extend the region in the 3' direction... This will be run streaming through a BAM file, so we should start
		// at the 5' end for all regions anyway.
		
		GenomeSpan from = new GenomeSpan(read.getReferenceName(), read.getAlignmentStart()-1, read.getAlignmentEnd()+extend);

		// the "to" is only an estimate of the proper span for now... when we parse the mate (later in the file), this will get fixed.
		// We want to extend the target region to make it more likely to merge other reads indicative of the SV 
    	GenomeSpan to = new GenomeSpan(read.getMateReferenceName(), read.getMateAlignmentStart()-1-extend, read.getMateAlignmentStart()+extend);
    
    	addTranslocation(read.getReadName(), from, to, read.getReadNegativeStrandFlag());
	}

	public void addTranslocation(String readName, GenomeSpan from, GenomeSpan to, boolean isNegStrand) {
		// We will only be scanning in the 5->3 direction, and the BAM file is sorted,
		// so we store a current translocation and use it to keep state. This also makes
		// processing significantly faster.
		//
		// If the new read doesn't fit the current, then we can assume it's part of a new
		// block.

		boolean found = false;
		List<Translocation> newbuf = new ArrayList<Translocation> ();
		
		for (Translocation current: buffer) {
			if (current.containsRead(readName)) {
				// Is this read already known?
				// (this shouldn't be common with the temp buffer stream approach)
				current.extendFrom(from);
				current.extendTo(to);
				current.addEvidence(readName, isNegStrand);
				newbuf.add(current);
				found = true;
			} else if (current.overlapsFrom(from)) {
				// do the "froms" overlap? 
				if (current.overlapsTo(to)) {
					// if the "tos" overlap, add this read
					current.extendFrom(from);
					current.extendTo(to);
					current.addEvidence(readName, isNegStrand);
					found = true;
				}
				// since we matched "froms", this region is still in play for the next read
				newbuf.add(current);
			} else {
				// if we don't overlap the "from", we never will... kick this one out or add it to the list.
				if (current.getEvidenceCount() >= minEvidence) {
					// only add if there is enough evidence (for this half of the SV, check the reciprocal later) 
					addAnnotation(current.getFrom(), current);
				}
			}
		}

		buffer = newbuf;

		if (!found) {
			Translocation current = new Translocation(from, to);
			current.addEvidence(readName, isNegStrand);
			buffer.add(current);
		}
	}
	
	public List<Translocation> peekBuffer() {
		return Collections.unmodifiableList(buffer);
	}
	
	public void done() {
		for (Translocation current: buffer) {
			if (current.getEvidenceCount() >= minEvidence) {
				addAnnotation(current.getFrom(), current);
			}
		}
		buffer.clear();
	}
	
	public void dump() {
		for (GenomeAnnotation<Translocation> t: IterUtils.wrap(annotations.iterator())) {
			if (t.getValue().getEvidenceCount() >= minEvidence) {
				if (t.getValue().getReciprocal() != null) {
					System.out.println(t.getValue() + "\t" + t.getValue().getReciprocal());
				}
			}
		}
	}

	public void findReciprocalMatches() {
//		Set<Translocation> matches = new HashSet<Translocation>();
		PrintStream out;
		try {
			out = new PrintStream(new FileOutputStream("tmp.txt"));
		} catch (FileNotFoundException e) {
			return;
		}
		
		for (GenomeAnnotation<Translocation> t: IterUtils.wrap(ProgressUtils.getIterator("Finding reciprocal matches", annotations.iterator(), new IncrementingStats(annotations.size())))) {
			if (t.getValue().getReciprocal() != null) {
				continue;
			}
			Translocation target = null;
			out.println("Trans t: " + t.getValue().getFrom()+ " => " + t.getValue().getTo());
			// look for all possible translocations in the "to" region
			for (Translocation potential : findAnnotation(t.getValue().getTo())) {
				out.print("  potential: " + potential.getFrom()+ " => " + potential.getTo());
				// for each of these, is there a reciprocal match?
				if (t.getValue().getFrom().overlaps(potential.getTo())) {
					out.print(" RECIP MATCH!");
					if (target == null) {
						out.println("");
						target = potential;
					} else {
						target.combine(potential);
						out.println(" Merging to: "+ target.getFrom()+ " => " + target.getTo());
					}
				} else {
					out.println(" NO RECIP MATCH!");
				}
			}
			if (target != null) {
				t.getValue().setReciprocal(target);
			}
		}
		out.close();
	}
}