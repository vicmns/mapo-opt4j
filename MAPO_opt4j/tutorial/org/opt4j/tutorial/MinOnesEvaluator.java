package org.opt4j.tutorial;

import org.opt4j.core.Objective;
import org.opt4j.core.Objective.Sign;
import org.opt4j.core.Objectives;
import org.opt4j.core.problem.Evaluator;

public class MinOnesEvaluator implements Evaluator<MinOnesResult> {

	Objective ones = new Objective("ones", Sign.MIN);

	public Objectives evaluate(MinOnesResult minOnesResult) {

		int value = 0;
		for (Boolean v : minOnesResult) {
			if (v!= null && v) {
				value++;
			}
		}

		Objectives objectives = new Objectives();
		objectives.add(ones, value);

		return objectives;
	}
}
