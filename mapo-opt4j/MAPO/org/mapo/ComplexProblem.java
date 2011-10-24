package org.mapo;

import org.opt4j.start.Constant;

import com.google.inject.Inject;

public class ComplexProblem {

	protected final int idxS;
	protected final int idxE;
	protected final String Gene; //Correct cDNA/DNA Sequence
	protected final String AB; //Mutated cDNA/RNA
	
	@Inject
	public ComplexProblem(@Constant(value = "srcGene", namespace = ComplexProblem.class) String srcGene, 
			@Constant(value = "mutadedGene", namespace = ComplexProblem.class) String mutadedGene, 
			@Constant(value = "idxS", namespace = ComplexProblem.class) int idxS,
			@Constant(value = "idxE", namespace = ComplexProblem.class) int idxE){
		this.Gene=srcGene;
		this.AB=mutadedGene;
		this.idxS=idxS;
		this.idxE=idxE;
	}
	
	public int getIdxS(){
		return this.idxS;
	}
	
	public int getIdxE(){
		return this.idxE;
	}
	
	public int getGeneLength(){
		return Gene.length();
	}
	
	public int getABLength(){
		return AB.length();
	}
	
	public String getGene(){
		return Gene;
	}
	
	public String getAB(){
		return AB;
	}
}
