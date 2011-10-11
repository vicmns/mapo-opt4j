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

package org.opt4j.optimizer.rs;

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
import org.opt4j.start.Constant;

import com.google.inject.Inject;

/**
 * The {@link RandomSearch} simply generates random {@link Individual}s and
 * evaluates them.
 * 
 * @author lukasiewycz
 * 
 */
public class RandomSearch extends AbstractOptimizer {

	protected final int batchsize;
	private final IndividualFactory individualFactory;

	/**
	 * Constructs a {@link RandomSearch}.
	 * 
	 * @param population
	 *            the population
	 * @param archive
	 *            the archive
	 * @param individualFactory
	 *            the individual creator
	 * @param completer
	 *            the completer
	 * @param control
	 *            the control
	 * @param iteration
	 *            the iteration counter
	 * @param batchsize
	 *            the size of the batch for an evaluation
	 */
	@Inject
	public RandomSearch(Population population, Archive archive, IndividualFactory individualFactory,
			IndividualCompleter completer, Control control, Iteration iteration,
			@Constant(value = "batchsize", namespace = RandomSearch.class) int batchsize) {
		super(population, archive, completer, control, iteration);
		this.batchsize = batchsize;
		this.individualFactory = individualFactory;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.core.optimizer.Optimizer#optimize()
	 */
	@Override
	public void optimize() throws StopException, TerminationException {
		while (iteration.value() < iteration.max()) {
			for (int i = 0; i < batchsize; i++) {
				Individual individual = individualFactory.create();
				population.add(individual);
			}
			nextIteration();
			population.clear();
		}
	}

}
