/**
 * Opt4J is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 * 
 * Opt4J is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with Opt4J. If not, see http://www.gnu.org/licenses/. 
 */
package org.opt4j.benchmark.knapsack;

import static java.lang.Math.max;
import static org.opt4j.core.Objective.INFEASIBLE;
import static org.opt4j.core.Objective.RANK_ERROR;
import static org.opt4j.core.Objective.Sign.MAX;
import static org.opt4j.core.Objective.Sign.MIN;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.inject.Inject;

import org.opt4j.core.Objective;
import org.opt4j.core.Objectives;
import org.opt4j.core.problem.Evaluator;

/**
 * The {@link KnapsackEvaluator} evaluates a given {@link ItemSelection} using
 * one {@link Objective} to sum up all knapsack overloads and one
 * {@link Objective} for the profit of each knapsack.
 * 
 * @author reimann, lukasiewycz
 * 
 */
public class KnapsackEvaluator implements Evaluator<ItemSelection> {
	protected Objective objectiveOverload = new Objective("overload", MIN, RANK_ERROR);
	protected final Map<Knapsack, Objective> profitObjectives = new HashMap<Knapsack, Objective>();
	protected KnapsackProblem problem;

	/**
	 * Creates a new {@link KnapsackEvaluator}.
	 * 
	 * @param problem
	 */
	@Inject
	public KnapsackEvaluator(KnapsackProblem problem) {
		this.problem = problem;
		int i = 0;
		for (Knapsack knapsack : problem.getKnapsacks()) {
			profitObjectives.put(knapsack, new Objective("profit" + i++, MAX));
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.core.problem.Evaluator#evaluate(org.opt4j.core.Phenotype)
	 */
	@Override
	public Objectives evaluate(ItemSelection itemSelection) {
		Objectives objectives = new Objectives();

		// evaluate overlead
		double overload = 0;
		for (Knapsack knapsack : problem.getKnapsacks()) {
			overload += getOverload(knapsack, itemSelection);
		}
		objectives.add(objectiveOverload, overload);
		objectives.setFeasible(overload == 0);
		
		// evaluate profit
		for (Knapsack knapsack : problem.getKnapsacks()) {
			if (objectives.isFeasible()) {
				objectives.add(profitObjectives.get(knapsack), getProfit(knapsack, itemSelection));
			} else {
				objectives.add(profitObjectives.get(knapsack), INFEASIBLE);
			}
		}
		return objectives;
	}

	/**
	 * Computes how much the given item exceed the defined capacity of the
	 * knapsack.
	 * 
	 * @param knapsack
	 *            the knapsack
	 * @param items
	 *            the items
	 * @return the overload
	 */
	protected double getOverload(Knapsack knapsack, Collection<Item> items) {
		double weight = 0;
		for (Item i : items) {
			weight += knapsack.getWeight(i);
		}
		return max(0, weight - knapsack.getCapacity());
	}

	/**
	 * Computes the profit for the given items for this knapsack.
	 * 
	 * @param knapsack
	 *            the knapsack
	 * @param items
	 *            the items
	 * @return the profit
	 */
	protected int getProfit(Knapsack knapsack, Collection<Item> items) {
		int profit = 0;
		for (Item i : items) {
			profit += knapsack.getProfit(i);
		}
		return profit;
	}
}
