package org.mapo;

import org.opt4j.core.problem.Decoder;
import org.opt4j.core.problem.PhenotypeWrapper;
import org.opt4j.genotype.StringGenotype;

public class ComplexDecoder implements Decoder<StringGenotype<String>, PhenotypeWrapper<String>> {
	@Override
	public PhenotypeWrapper<String> decode(StringGenotype<String> genotype) {
		String s="";
		s=genotype.toString();
		return new PhenotypeWrapper<String>(s);
	}
}
