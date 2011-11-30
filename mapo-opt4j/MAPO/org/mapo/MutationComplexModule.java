package org.mapo;

import org.opt4j.config.annotations.Info;
import org.opt4j.config.annotations.Required;
import org.opt4j.operator.mutate.AdaptiveMutationRate;
import org.opt4j.operator.mutate.ConstantMutationRate;
import org.opt4j.operator.mutate.MutateComplexInversion;
import org.opt4j.operator.mutate.MutateModule;
import org.opt4j.operator.mutate.MutationRate;
import org.opt4j.operator.mutate.BasicMutateModule.MutationRateType;
import org.opt4j.start.Constant;

@Info("Setting for the complex mutate classOperators for genotype variation.")
public class MutationComplexModule extends MutateModule{
	
	@Info("The type of mutation rate.")
	protected MutationRateType mutationRateType = MutationRateType.CONSTANT;
	
	@Required(property = "mutationRateType", elements = { "CONSTANT" })
	@Constant(value = "rate", namespace = ConstantMutationRate.class)
	protected double mutationRate = 0.01;
	
	@Info("The type of mutation to use in a DX/TX-COMPLEX.")
	protected ComplexMutation complexMutation = ComplexMutation.INVERSION;
	
	/**
	 * The type of mutation to use in a DX/TX-COMPLEX
	 * @author vicmns
	 *
	 */
	public enum ComplexMutation{
		/**
		 * Select two random and inverts the values.
		 */
		INVERSION,
		/**
		 * Select a random close structure, then rest a random
		 * value for any non-close structure adjacent.
		 */
		CLOSECOMPLEMENT,
		/**
		 * Select a random non close structure, the select an adjacent close
		 * structure and delete both.
		 */
		DELETION;
	}
	
	/**
	 * The type of the used mutation rate.
	 * 
	 * @author lukasiewycz
	 * 
	 */
	public enum MutationRateType {
		/**
		 * Use a constant mutation rate.
		 */
		CONSTANT,
		/**
		 * Use a variable (adaptive) mutation rate.
		 */
		ADAPTIVE;
	}
	
	public ComplexMutation getComplexMutation(){
		return complexMutation;
	}
	
	public void setComplexMutation(ComplexMutation complexMutation){
		this.complexMutation=complexMutation;
	}
	
	/**
	 * Returns the {@link MutationRateType}.
	 * 
	 * @return the type of mutation rate
	 */
	public MutationRateType getMutationRateType() {
		return mutationRateType;
	}

	/**
	 * Sets the {@link MutationRateType}.
	 * 
	 * @param mutationRateType
	 *            the type of mutation rate
	 */
	public void setMutationRateType(MutationRateType mutationRateType) {
		this.mutationRateType = mutationRateType;
	}
	
	/**
	 * Returns the {@link MutationRate} as a double value.
	 * 
	 * @return the mutation rate
	 */
	public double getMutationRate() {
		return mutationRate;
	}

	/**
	 * Sets the {@link MutationRate}.
	 * 
	 * @param mutationRate
	 *            the mutation rate
	 */
	public void setMutationRate(double mutationRate) {
		this.mutationRate = mutationRate;
	}
	
	@Override
	protected void config(){
		switch (mutationRateType) {
		case ADAPTIVE:
			bind(MutationRate.class).to(AdaptiveMutationRate.class).in(SINGLETON);
			break;
		case CONSTANT:
			bind(MutationRate.class).to(ConstantMutationRate.class).in(SINGLETON);
		}
		
		switch(complexMutation){
		case INVERSION:
			//bind(MutateComplex.class).to(MutateComplexInversion.class).in(SINGLETON);
			addOperator(MutateComplexInversion.class);
			break;
		case CLOSECOMPLEMENT:
			
			break;
		case DELETION:
			
			break;
		default:
			addOperator(MutateComplexInversion.class);
			break;
		}
	}

}
