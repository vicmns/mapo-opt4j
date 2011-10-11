package org.opt4j.tutorial;

import java.util.Collection;

import org.opt4j.core.Genotype;
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
import org.opt4j.operator.copy.Copy;
import org.opt4j.operator.mutate.Mutate;
import org.opt4j.optimizer.ea.Selector;

import com.google.inject.Inject;

public class MutateOptimizer extends AbstractOptimizer {

	protected final IndividualFactory individualFactory;
	
	protected final Mutate<Genotype> mutate;

	protected final Copy<Genotype> copy;

	protected final Selector selector;

	public static final int POPSIZE = 100;

	public static final int OFFSIZE = 25;
	

	@Inject
	public MutateOptimizer(Population population, Archive archive, IndividualFactory individualFactory,
			IndividualCompleter completer, Control control, Selector selector, Mutate<Genotype> mutate,
			Copy<Genotype> copy, Iteration iteration) {
		super(population, archive, completer, control, iteration);
		this.individualFactory = individualFactory;
		this.mutate = mutate;
		this.copy = copy;
		this.selector = selector;
	}

	public void optimize() throws TerminationException, StopException {
		selector.init(OFFSIZE + POPSIZE);

		for (int i = 0; i < 100; i++) {
			population.add(individualFactory.create());
		}

		nextIteration();

		while (iteration.value() < iteration.max()) {

			Collection<Individual> parents = selector.getParents(OFFSIZE, population);

			for (Individual parent : parents) {
				Genotype genotype = copy.copy(parent.getGenotype());
				mutate.mutate(genotype, 0.1);

				Individual child = individualFactory.create(genotype);
				population.add(child);
			}

			completer.complete(population);

			Collection<Individual> lames = selector.getLames(OFFSIZE, population);
			population.removeAll(lames);

			nextIteration();

		}

	}

}
