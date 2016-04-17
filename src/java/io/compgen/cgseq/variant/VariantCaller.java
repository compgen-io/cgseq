package io.compgen.cgseq.variant;

import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCall;

import java.util.List;
import java.util.Map;

public interface VariantCaller {
	public abstract VariantResults calcVariant(List<PileupBaseCall> calls, String ref);
	public abstract Map<String, String> getInfoFields();
	public abstract Map<String, String> getFormatFields();

}