package org.mapo;

import java.util.ArrayList;

import org.opt4j.core.Phenotype;
import org.opt4j.core.problem.Decoder;
import org.opt4j.genotype.DynamicListGenotype;

public class ComplexDecoder implements Decoder<DynamicListGenotype<Object>, ComplexPhenotype> {
	
	@Override
	public ComplexPhenotype decode(DynamicListGenotype<Object> genotype) {
		ComplexPhenotype complex = new ComplexPhenotype();
		complex.addAll(genotype);
		return complex;
	}
}
