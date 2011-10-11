package org.mapo;

import org.opt4j.core.Phenotype;
import org.opt4j.core.problem.Decoder;
import org.opt4j.core.problem.PhenotypeWrapper;
import org.opt4j.genotype.SelectGenotype;

public class complexDecoder implements Decoder<SelectGenotype<String>, PhenotypeWrapper<String>> {
	@Override
	public PhenotypeWrapper<String> decode(SelectGenotype<String> genotype) {
		String s="";
		for (int i = 0; i < genotype.size(); i++) {
			s += genotype.getValue(i);
		}
		return new PhenotypeWrapper<String>(s);
	}
}
