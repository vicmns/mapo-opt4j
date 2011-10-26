package org.mapo;

import org.opt4j.operator.mutate.MutateComplexInversion;
import org.opt4j.operator.mutate.MutateModule;

public class MutationComplexModule extends MutateModule{
	
	@Override
	protected void config(){
		addOperator(MutateComplexInversion.class);
	}

}
