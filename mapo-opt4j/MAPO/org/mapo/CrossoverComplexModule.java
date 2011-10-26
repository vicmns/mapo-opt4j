package org.mapo;

import org.opt4j.operator.crossover.CrossoverComplexOnMutation;
import org.opt4j.operator.crossover.CrossoverModule;

public class CrossoverComplexModule extends CrossoverModule{
	
	@Override
	protected void config(){
		addOperator(CrossoverComplexOnMutation.class);
	}
}
