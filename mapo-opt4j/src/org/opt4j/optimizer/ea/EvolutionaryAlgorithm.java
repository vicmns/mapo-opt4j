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
package org.opt4j.optimizer.ea;

import java.util.Collection;

import org.opt4j.core.Individual;
import org.opt4j.core.IndividualFactory;
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
import org.opt4j.start.Constant;

import com.google.inject.Inject;

/**
 * The {@link EvolutionaryAlgorithm} is an implementation of an Evolutionary
 * Algorithm based on the operators {@link Crossover} and {@link Mutate}. It
 * uses a {@link Selector} for the mating and environmental selection.
 * 
 * @author lukasiewycz, glass
 * 
 */
public class EvolutionaryAlgorithm extends AbstractOptimizer {

	protected final int alpha;

	protected final int lambda;

	protected final int mu;

	protected final Selector selector;

	protected final Mating mating;

	private final IndividualFactory individualFactory;

	/**
	 * Constructs an {@link EvolutionaryAlgorithm} with a {@link Population}, an
	 * {@link Archive}, an {@link IndividualFactory}, a
	 * {@link IndividualCompleter}, a {@link Control}, a {@link Selector}, a
	 * {@link Mating}, the number of generations, the population size, the
	 * number of parents, the number of offspring, and a random number
	 * generator.
	 * 
	 * @param population
	 *            the population
	 * @param archive
	 *            the archive
	 * @param individualFactory
	 *            the individual factory
	 * @param completer
	 *            the completer
	 * @param control
	 *            the control
	 * @param selector
	 *            the selector
	 * @param mating
	 *            the mating
	 * @param iteration
	 *            the iteration counter
	 * @param alpha
	 *            the population size
	 * @param mu
	 *            the number of parents
	 * @param lambda
	 *            the number of offspring
	 */
	@Inject
	public EvolutionaryAlgorithm(Population population, Archive archive, IndividualFactory individualFactory,
			IndividualCompleter completer, Control control, Selector selector, Mating mating, Iteration iteration,
			@Constant(value = "alpha", namespace = EvolutionaryAlgorithm.class) int alpha,
			@Constant(value = "mu", namespace = EvolutionaryAlgorithm.class) int mu,
			@Constant(value = "lambda", namespace = EvolutionaryAlgorithm.class) int lambda) {
		super(population, archive, completer, control, iteration);
		this.selector = selector;
		this.mating = mating;
		this.alpha = alpha;
		this.mu = mu;
		this.lambda = lambda;
		this.individualFactory = individualFactory;

		if (alpha <= 0) {
			throw new IllegalArgumentException("Invalid alpha: " + alpha);
		}

		if (mu <= 0) {
			throw new IllegalArgumentException("Invalid mu: " + mu);
		}

		if (lambda <= 0) {
			throw new IllegalArgumentException("Invalid lambda: " + lambda);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.core.optimizer.Optimizer#optimize()
	 */
	@Override
	public void optimize() throws TerminationException, StopException {

		selector.init(alpha + lambda);

		while (iteration.value() < iteration.max()) {

			// set the number of new offspring individuals per generation
			int offspringCount = lambda;

			// 1) fill the population until it reaches the size alpha
			while (population.size() < alpha && offspringCount > 0) {
				population.add(individualFactory.create());
				offspringCount--;
			}

			// 2) generate offspring by mating
			if (offspringCount > 0) {
				if (offspringCount < lambda) { // evaluate new individuals first
					completer.complete(population);
				}

				Collection<Individual> parents = selector.getParents(mu, population);
				Collection<Individual> offspring = mating.getOffspring(offspringCount, parents);
				population.addAll(offspring);
			}

			// evaluate offspring before selecting lames
			completer.complete(population);

			if (population.size() > alpha) {
				Collection<Individual> lames = selector.getLames(population.size() - alpha, population);
				population.removeAll(lames);
			}

			nextIteration();
		}
	}
}
