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

package org.opt4j.optimizer.sa;

import java.util.Random;

import org.opt4j.common.random.Rand;
import org.opt4j.core.Genotype;
import org.opt4j.core.Individual;
import org.opt4j.core.IndividualFactory;
import org.opt4j.core.Objectives;
import org.opt4j.core.optimizer.AbstractOptimizer;
import org.opt4j.core.optimizer.Archive;
import org.opt4j.core.optimizer.Control;
import org.opt4j.core.optimizer.IndividualCompleter;
import org.opt4j.core.optimizer.Iteration;
import org.opt4j.core.optimizer.Population;
import org.opt4j.core.optimizer.StopException;
import org.opt4j.core.optimizer.TerminationException;
import org.opt4j.operator.copy.Copy;
import org.opt4j.operator.neighbor.Neighbor;

import com.google.inject.Inject;

/**
 * The {@link SimulatedAnnealing} is a standard implementation of the
 * optimization heuristic. This heuristic optimizes a single objective
 * (multi-objective problems are optimized by the sum of all objectives).
 * 
 * @author lukasiewycz
 * 
 */
public class SimulatedAnnealing extends AbstractOptimizer {

	protected final Random random;

	protected final Neighbor<Genotype> neighbor;

	protected final Copy<Genotype> copy;

	protected final CoolingSchedule coolingSchedule;

	private final IndividualFactory individualFactory;

	/**
	 * Constructs a new {@link SimulatedAnnealing}.
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
	 * @param random
	 *            the random number generator
	 * @param neighbor
	 *            the neighbor operator
	 * @param copy
	 *            the copy operator
	 * @param iteration
	 *            the iteration counter
	 * @param coolingSchedule
	 *            the cooling schedule
	 */
	@Inject
	public SimulatedAnnealing(Population population, Archive archive, IndividualFactory individualFactory,
			IndividualCompleter completer, Control control, Rand random, Neighbor<Genotype> neighbor,
			Copy<Genotype> copy, Iteration iteration, CoolingSchedule coolingSchedule) {
		super(population, archive, completer, control, iteration);
		this.random = random;
		this.neighbor = neighbor;
		this.copy = copy;
		this.coolingSchedule = coolingSchedule;
		this.individualFactory = individualFactory;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.core.optimizer.Optimizer#optimize()
	 */
	@Override
	public void optimize() throws TerminationException, StopException {

		Individual x = individualFactory.create();
		Individual y = null;

		population.add(x);
		completer.complete(population);
		archive.update(x);
		nextIteration();

		double fx = f(x);
		double fy;

		while (iteration.value() < iteration.max()) {

			Genotype g = copy.copy(x.getGenotype());
			neighbor.neighbor(g);

			y = individualFactory.create(g);

			completer.complete(y);
			archive.update(y);

			fy = f(y);
			// boolean value that indicates a switch of the individuals
			boolean sw = false;

			if (fy <= fx) {
				sw = true;
			} else {
				double a = (fx - fy) / coolingSchedule.getTemperature(iteration.value(), iteration.max());
				double e = Math.exp(a);
				if (random.nextDouble() < e) {
					sw = true;
				}
			}

			if (sw) {
				population.remove(x);
				population.add(y);
				fx = fy;
				x = y;
			}

			nextIteration();
		}
	}

	/**
	 * Calculates the sum of the {@link Objectives} of one {@link Individual}.
	 * 
	 * @param individual
	 *            the individual
	 * @return the sum of the objective values
	 */
	protected double f(Individual individual) {
		Objectives objectives = individual.getObjectives();

		double sum = 0;

		for (double d : objectives.array()) {
			sum += d;
		}

		return sum;
	}

}
