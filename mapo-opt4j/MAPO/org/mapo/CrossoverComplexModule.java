package org.mapo;

import org.opt4j.config.annotations.Info;
import org.opt4j.operator.crossover.CrossoverComplexOnMutation;
import org.opt4j.operator.crossover.CrossoverComplexUniform;
import org.opt4j.operator.crossover.CrossoverModule;
import org.opt4j.operator.crossover.BasicCrossoverModule.DynamicListType;

@Info("Setting for the complex crossover classOperators for genotype variation.")
public class CrossoverComplexModule extends CrossoverModule{
	
	@Info("The type of the crossover operator for the Complex genotype.")
	protected ComplexCrossover complexCrossover = ComplexCrossover.SBSCROSSOVER;
	
	public enum ComplexCrossover{
		SBSCROSSOVER,
		UNIFORMCROSSOVER;
	}
	
	public ComplexCrossover getComplexCrossover(){
		return complexCrossover;
	}
	public void setComplexCrossover(ComplexCrossover complexCrossover){
		this.complexCrossover=complexCrossover;
	}
	
	@Override
	protected void config(){
		switch(complexCrossover){
			case SBSCROSSOVER:
				addOperator(CrossoverComplexOnMutation.class);
				break;
			case UNIFORMCROSSOVER:
				addOperator(CrossoverComplexUniform.class);
				break;
			default:
				addOperator(CrossoverComplexOnMutation.class);
				break;
		}
	}
}
