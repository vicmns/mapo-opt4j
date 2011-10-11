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

package org.opt4j.common.archive;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.opt4j.core.Individual;
import org.opt4j.optimizer.ea.Nsga2;
import org.opt4j.start.Constant;

import com.google.inject.Inject;

/**
 * The {@link CrowdingArchive} is based on the crowding distance of the
 * {@link Nsga2}.
 * 
 * @author lukasiewycz
 * 
 */
public class CrowdingArchive extends BoundedArchive {

	/**
	 * Constructs a {@link CrowdingArchive}.
	 * 
	 * @param capacity
	 *            capacity of this archive (using namespace
	 *            {@link BoundedArchive})
	 */
	@Inject
	public CrowdingArchive(@Constant(value = "capacity", namespace = BoundedArchive.class) int capacity) {
		super(capacity);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.opt4j.common.archive.AbstractArchive#updateWithNondominated(java.
	 * util.Collection)
	 */
	@Override
	protected boolean updateWithNondominated(Collection<Individual> candidates) {
		boolean changed = false;
		if (this.size() + candidates.size() <= capacity) {
			changed = addCheckedIndividuals(candidates);
		} else {
			candidates.addAll(this);
			Map<Individual, Double> map = getValues(candidates);
			List<Individual> list = order(map);
			List<Individual> worst = list.subList(capacity, list.size());
			candidates.removeAll(worst);

			this.retainAll(candidates);
			for (Individual i : candidates) {
				if (!this.contains(i)) {
					changed |= addCheckedIndividual(i);
				}
			}
		}
		return changed;
	}

	/**
	 * Calculates the crowding distance for a collection of {@link Individual}s.
	 * 
	 * @param individuals
	 *            the individuals
	 * @return a map of the individuals to the corresponding crowding distance
	 */
	public Map<Individual, Double> getValues(Collection<Individual> individuals) {
		Map<Individual, Double> map = new LinkedHashMap<Individual, Double>();
		List<Individual> list = new ArrayList<Individual>(individuals);

		for (Individual individual : list) {
			map.put(individual, 0.0);
		}

		int m = list.get(0).getObjectives().array().length;

		for (int i = 0; i < m; i++) {
			final int d = i;
			Collections.sort(list, new Comparator<Individual>() {
				@Override
				public int compare(Individual o1, Individual o2) {
					double i1 = o1.getObjectives().array()[d];
					double i2 = o2.getObjectives().array()[d];
					return (int) Math.signum(i1 - i2);
				}
			});

			map.put(list.get(0), Double.MAX_VALUE);
			map.put(list.get(list.size() - 1), Double.MAX_VALUE);

			for (int j = 1; j < list.size() - 1; j++) {
				double p = list.get(j - 1).getObjectives().array()[i];
				double n = list.get(j + 1).getObjectives().array()[i];
				Individual ind = list.get(j);
				map.put(ind, map.get(ind) + (n - p));
			}
		}

		return map;
	}

	/**
	 * Returns an ordered list of the {@link Individual}s corresponding to their
	 * crowding distance. {@link Individual}s with a high crowding distance are
	 * sorted to the front.
	 * 
	 * @param values
	 *            the map of individuals to their crowding distance values
	 * @return the sorted list of individuals based on their crowding distance
	 */
	public List<Individual> order(final Map<Individual, Double> values) {
		List<Individual> list = new ArrayList<Individual>(values.keySet());
		Collections.sort(list, new Comparator<Individual>() {
			@Override
			public int compare(Individual o1, Individual o2) {
				double i1 = values.get(o1);
				double i2 = values.get(o2);
				return (int) Math.signum(i2 - i1);
			}
		});
		return list;
	}

}
