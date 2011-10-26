package org.opt4j.operator.crossover;

import java.util.Random;

import org.opt4j.genotype.DynamicListGenotype;

import com.google.inject.Inject;

public class CrossoverComplexOnMutation extends CrossoverComplexOnePoint<DynamicListGenotype<?>> implements CrossoverComplex {
	
	@Inject
	public CrossoverComplexOnMutation(Random random){
		super(random);
	}
}
