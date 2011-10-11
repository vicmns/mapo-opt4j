package org.opt4j.optimizer.ea;

import static org.opt4j.core.Objective.Sign.MAX;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.opt4j.core.Individual;
import org.opt4j.core.IndividualStateListener;
import org.opt4j.core.Objective;
import org.opt4j.core.Objectives;
import org.opt4j.core.Value;

import com.google.inject.Singleton;

@Singleton
public class Normalizer implements IndividualStateListener {
	private final Map<Objective, Double> minValues = new HashMap<Objective, Double>();
	private final Map<Objective, Double> maxValues = new HashMap<Objective, Double>();

	/**
	 * Returns normalized objectives. Each objective is in the range between 0
	 * and 1 and have to be minimized. Here, 0 is the smallest value seen so far
	 * for this objective for all evaluated {@link Individual}s and 1 the
	 * biggest value, respectively. If an objective is infeasible, it is set to
	 * 1.
	 * 
	 * @param objectives
	 *            the objectives to normalize
	 * @return the normalized objectives
	 */
	public Objectives normalize(Objectives objectives) {
		Objectives normalized = new Objectives();

		for (Entry<Objective, Value<?>> entry : objectives) {
			Objective objective = entry.getKey();
			double oldvalue = toMinProblem(entry);

			double newvalue = 1;

			if (oldvalue != Double.MAX_VALUE) {
				double min = minValues.get(objective);
				double max = maxValues.get(objective);
				newvalue = (oldvalue - min) / (max - min);
			}
			normalized.add(objective, newvalue);
		}

		return normalized;
	}

	@Override
	public void inidividualStateChanged(Individual individual) {
		if (individual.isEvaluated()) {
			for (Entry<Objective, Value<?>> entry : individual.getObjectives()) {
				Objective objective = entry.getKey();
				double value = toMinProblem(entry);
				if (value != Double.MAX_VALUE) {
					if (minValues.get(objective) == null || value < minValues.get(objective)) {
						minValues.put(objective, value);
					}
					if (maxValues.get(objective) == null || value > maxValues.get(objective)) {
						maxValues.put(objective, value);
					}
				}
			}
		}
	}

	public static final double toMinProblem(Entry<Objective, Value<?>> entry) {
		Objective objective = entry.getKey();
		Value<?> value = entry.getValue();

		Double v = value.getDouble();

		Double result = null;
		if (v == null) {
			result = Double.MAX_VALUE;
		} else if (objective.getSign() == MAX) {
			result = -v;
		} else {
			result = v;
		}

		assert result != null;
		return result;
	}
}
