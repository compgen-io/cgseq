package io.compgen.cgseq.variant;

import io.compgen.ngsutils.pileup.PileupRecord.PileupBaseCall;

import java.util.List;

public interface VariantCaller {
	public abstract VariantResults calcVariant(List<PileupBaseCall> calls, String ref);
	public abstract List<String> getInfoFields();
	public abstract List<String> getFormatFields();
	public abstract String getInfoFieldDescription(String k);
	public abstract String getFormatFieldDescription(String k);

}