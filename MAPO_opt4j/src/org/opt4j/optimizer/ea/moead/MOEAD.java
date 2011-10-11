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
package org.opt4j.optimizer.ea.moead;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.Vector;

import org.opt4j.common.random.Rand;
import org.opt4j.core.Individual;
import org.opt4j.core.IndividualFactory;
import org.opt4j.core.Objective;
import org.opt4j.core.Objectives;
import org.opt4j.core.Value;
import org.opt4j.core.optimizer.AbstractOptimizer;
import org.opt4j.core.optimizer.Archive;
import org.opt4j.core.optimizer.Control;
import org.opt4j.core.optimizer.IndividualCompleter;
import org.opt4j.core.optimizer.Iteration;
import org.opt4j.core.optimizer.Population;
import org.opt4j.core.optimizer.StopException;
import org.opt4j.core.optimizer.TerminationException;
import org.opt4j.operator.crossover.Crossover;
import org.opt4j.operator.mutate.Mutate;
import org.opt4j.optimizer.ea.Mating;
import org.opt4j.optimizer.ea.Normalizer;
import org.opt4j.optimizer.ea.Selector;
import org.opt4j.start.Constant;

import com.google.inject.Inject;
import com.google.inject.Provider;

/**
 * The {@link MOEAD} is an implementation of an Evolutionary Algorithm based on
 * the operators {@link Crossover} and {@link Mutate}. It uses decomposition for
 * the mating and environmental selection.
 * 
 * @see "Q. Zhang: MOEA/D: A Multiobjective Evolutionary Algorithm Based on
 *      Decomposition. IEEE Transactions on Evolutionary Computation, vol. 11,
 *      no. 6, pp. 712-731, Dec. 2007."
 * 
 * @author reimann
 * 
 */
public class MOEAD extends AbstractOptimizer {

	private final int subproblemsN;

	private final int neighborsT;

	private final Mating mating;

	private final Rand random;

	private final Normalizer normalizer;

	private final Provider<Subproblem> subproblemProvider;

	/**
	 * Constructs a {@link MOEAD} with a {@link Population}, an {@link Archive},
	 * an {@link IndividualFactory}, a {@link IndividualCompleter}, a
	 * {@link Control}, a {@link Selector}, a {@link Mating}, the number of
	 * generations, the population size, the number of parents, the number of
	 * offspring, and a random number generator.
	 * 
	 * @param population
	 *            the population
	 * @param archive
	 *            the archive
	 * @param completer
	 *            the completer
	 * @param control
	 *            the control
	 * @param mating
	 *            the mating
	 * @param iteration
	 *            the iteration counter
	 * @param subproblemsN
	 *            the population size
	 * @param neighborsT
	 *            the number of parents
	 */
	@Inject
	public MOEAD(Population population, Archive archive, IndividualCompleter completer, Control control, Mating mating,
			Iteration iteration, Rand random, Normalizer normalizer, Provider<Subproblem> subproblemProvider,
			@Constant(value = "n", namespace = MOEAD.class) int subproblemsN,
			@Constant(value = "t", namespace = MOEAD.class) int neighborsT) {
		super(population, archive, completer, control, iteration);
		this.subproblemProvider = subproblemProvider;
		this.mating = mating;
		this.subproblemsN = subproblemsN;
		this.neighborsT = neighborsT;
		this.random = random;
		this.normalizer = normalizer;

		if (subproblemsN <= 0) {
			throw new IllegalArgumentException("Invalid alpha: " + subproblemsN);
		}

		if (neighborsT <= 0) {
			throw new IllegalArgumentException("Invalid mu: " + neighborsT);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.core.optimizer.Optimizer#optimize()
	 */
	@Override
	public void optimize() throws TerminationException, StopException {

		List<Subproblem> subproblems = new ArrayList<Subproblem>();
		for (int i = 0; i < subproblemsN; i++) {
			Subproblem subproblem = subproblemProvider.get();
			subproblems.add(subproblem);
			population.add(subproblem.getBest());
		}

		completer.complete(population);

		for (Subproblem subproblem : subproblems) {
			subproblem.initWeightVector();
		}
		// List<Vector<Double>> vectors = initVectors();
		computeNeighborhood(subproblems);

		while (iteration.value() < iteration.max()) {

			for (Subproblem subproblem : subproblems) {
				// Reproduction
				List<Subproblem> neighbors = subproblem.getNeighbors();
				int k = random.nextInt(neighbors.size());
				int l = random.nextInt(neighbors.size());
				Individual parent1 = neighbors.get(k).getBest();
				Individual parent2 = neighbors.get(l).getBest();
				Individual y = mating.getOffspring(1, parent1, parent2).iterator().next();
				population.add(y);
				subproblem.setChild(y);
			}

			// evaluate offspring before selecting lames
			completer.complete(population);

			// updateReferencePoints(minValues, maxValues, ys);
			for (Subproblem subproblem : subproblems) {
				Individual individual = subproblem.getChild();
				for (Subproblem neighbor : subproblem.getNeighbors()) {
					double distanceY = distance(individual, neighbor.getScalar());
					double distanceX = distance(neighbor.getBest(), neighbor.getScalar());
					if (distanceY <= distanceX) {
						neighbor.setBest(individual);
					}
				}
			}
			nextIteration();
			population.clear();
			for (Subproblem subproblem : subproblems) {
				population.add(subproblem.getBest());
			}
		}
	}

	private double distance(Individual individual, Vector<Double> vector) {
		assert individual.isEvaluated();
		// normalize according to Z07 section V.G.2)
		Objectives objectives = normalizer.normalize(individual.getObjectives());

		int i = 0;
		double d = 0.0;
		for (Entry<Objective, Value<?>> o : objectives) {
			double value1 = o.getValue().getDouble();
			double value2 = vector.get(i);
			d = Math.max(value1 * value2, d);
			i++;
		}
		return d;
	}

	private void computeNeighborhood(List<Subproblem> subproblems) {
		for (Subproblem subproblem : subproblems) {
			TreeSet<Subproblem> orderedSet = new TreeSet<Subproblem>(new EuclideanDistance(subproblem));
			for (Subproblem otherproblem : subproblems) {
				orderedSet.add(otherproblem);
			}
			int j = 0;
			for (Subproblem neighbor : orderedSet) {
				if (j >= neighborsT) {
					break;
				}
				subproblem.addNeighbor(neighbor);
				j++;
			}
		}
	}

	/**
	 * The {@link EuclideanDistance} is a {@link Comparator} according to Zhang.
	 * 
	 * @author reimann
	 * 
	 */
	private class EuclideanDistance implements Comparator<Subproblem> {
		// TODO more efficient: compute distance o--root once
		private final Vector<Double> root;

		public EuclideanDistance(Subproblem rootProblem) {
			this.root = rootProblem.getScalar();
		}

		@Override
		public int compare(Subproblem subproblem1, Subproblem subproblem2) {
			Vector<Double> o1 = subproblem1.getScalar();
			Vector<Double> o2 = subproblem2.getScalar();
			assert o1.size() == o2.size();
			assert o1.size() == root.size();

			Double distance1 = 0.0;
			Double distance2 = 0.0;
			for (int i = 0; i < o1.size(); i++) {
				double dist1 = root.get(i) - o1.get(i);
				distance1 += dist1 * dist1;

				double dist2 = root.get(i) - o2.get(i);
				distance2 += dist2 * dist2;
			}

			return distance1.compareTo(distance2);
		}
	}

}
