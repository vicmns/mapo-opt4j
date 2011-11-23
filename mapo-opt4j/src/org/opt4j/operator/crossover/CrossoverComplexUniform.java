package org.opt4j.operator.crossover;

import java.util.Random;

import org.mapo.ComplexProblem;
import org.opt4j.genotype.DynamicListGenotype;

import com.google.inject.Inject;

public class CrossoverComplexUniform extends CrossoverComplexNPoint<DynamicListGenotype<?>> implements CrossoverComplex {
	
	@Inject
	public CrossoverComplexUniform(Random random, ComplexProblem problem){
		super(random, problem);
	}
	
}
