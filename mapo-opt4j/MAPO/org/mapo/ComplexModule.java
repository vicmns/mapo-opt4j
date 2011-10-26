package org.mapo;

import org.opt4j.core.problem.ProblemModule;

public class ComplexModule extends ProblemModule {
	public ComplexModule(){
		super();
	}
	
	public void config(){
		bind(ComplexProblem.class).in(SINGLETON);
		bindProblem(ComplexCreator.class, ComplexDecoder.class, ComplexEvaluator.class);
	}
}
