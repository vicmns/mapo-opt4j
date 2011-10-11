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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.opt4j.common.random.Rand;
import org.opt4j.core.Individual;
import org.opt4j.core.Objectives;
import org.opt4j.start.Constant;

import com.google.inject.Inject;

/**
 * The {@link Nsga2} {@link Selector}.
 * 
 * @see Nsga2Module
 * @see "A Fast Elitist Non-Dominated Sorting Genetic Algorithm for
 *      Multi-Objective Optimization: NSGA-II, K. Deb, Samir Agrawal, Amrit
 *      Pratap, and T. Meyarivan, Parallel MockProblem Solving from Nature,
 *      2000"
 * @author lukasiewycz, noorshams
 * 
 */
public class Nsga2 implements Selector {

	protected final Random random;
	protected final int tournament;

	protected final Map<Individual, Integer> map = new LinkedHashMap<Individual, Integer>();
	protected Individual[] ind = new Individual[0];
	protected int[] rank = new int[0];
	protected double[] dist = new double[0];

	protected Integer m = null;

	protected List<List<Integer>> fronts;

	/**
	 * Constructs a {@link Nsga2} {@link Selector}.
	 * 
	 * @param random
	 *            the random number generator
	 * @param tournament
	 *            the tournament value
	 */
	@Inject
	public Nsga2(Rand random, @Constant(value = "tournament", namespace = Nsga2.class) int tournament) {
		this.random = random;
		this.tournament = tournament;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.optimizer.ea.Selector#init(int)
	 */
	@Override
	public void init(int maxsize) {
		map.clear();

		ind = new Individual[maxsize];
		rank = new int[maxsize];
		dist = new double[maxsize];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.optimizer.ea.Selector#getParents(int,
	 * java.util.Collection)
	 */
	@Override
	public Collection<Individual> getParents(int mu, Collection<Individual> population) {
		if (synchronize(population)) {
			fronts();
		}
		List<Integer> all = new ArrayList<Integer>(map.values());
		List<Individual> parents = new ArrayList<Individual>();

		int size = all.size();

		// To avoid multiple distance calculations, save the calculated ranks
		// (front numbers)
		Set<Integer> alreadyCalculatedRanks = new LinkedHashSet<Integer>();

		for (int i = 0; i < mu; i++) {
			int winner = all.get(random.nextInt(size));

			for (int t = 0; t < tournament; t++) {
				int opponent = all.get(random.nextInt(size));
				if (rank[opponent] < rank[winner] || opponent == winner) {
					winner = opponent;
				} else if (rank[opponent] == rank[winner]) {
					// The winner is determined considering the crowding
					// distance

					List<Integer> front = fronts.get(rank[winner]);
					if (!(alreadyCalculatedRanks.contains(rank[winner]))) {
						// calculation needed
						alreadyCalculatedRanks.add(rank[winner]);
						calcDistance(front);
					}

					// Opponent wins, if it has a better crowding distance
					if (dist[opponent] > dist[winner]) {
						winner = opponent;
					}

				}
			}

			parents.add(ind[winner]);
		}

		return parents;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.optimizer.ea.Selector#getLames(int, java.util.Collection)
	 */
	@Override
	public Collection<Individual> getLames(int n, Collection<Individual> population) {
		synchronize(population);

		List<Integer> remove = new ArrayList<Integer>();
		int size = n;

		List<List<Integer>> fronts = fronts();
		Collections.reverse(fronts);

		for (List<Integer> front : fronts) {
			if (remove.size() + front.size() < size) {
				remove.addAll(front);
			} else {
				calcDistance(front);
				Collections.sort(front, new DistanceSort());
				remove.addAll(front.subList(0, size - remove.size()));
			}
		}

		List<Individual> truncate = new ArrayList<Individual>();
		for (Integer k : remove) {
			truncate.add(ind[k]);
		}

		return truncate;
	}

	/**
	 * Evaluate the fronts and set the correspondent rank values.
	 * 
	 * @return the fronts
	 */
	public List<List<Integer>> fronts() {

		List<Integer> pop = new ArrayList<Integer>(map.values());

		List<List<Integer>> fronts = new ArrayList<List<Integer>>();

		@SuppressWarnings("unchecked")
		List<Integer>[] S = new List[ind.length];
		int[] n = new int[ind.length];

		for (int e : pop) {
			S[e] = new ArrayList<Integer>();
			n[e] = 0;
		}

		for (int i = 0; i < pop.size(); i++) {
			for (int j = i + 1; j < pop.size(); j++) {
				int p = pop.get(i);
				int q = pop.get(j);

				Objectives po = ind[p].getObjectives();
				Objectives qo = ind[q].getObjectives();

				if (po.dominates(qo)) {
					S[p].add(q);
					n[q]++;
				} else if (qo.dominates(po)) {
					S[q].add(p);
					n[p]++;
				}
			}
		}

		List<Integer> f1 = new ArrayList<Integer>();
		for (int i : pop) {
			if (n[i] == 0) {
				f1.add(i);
			}
		}
		fronts.add(f1);
		List<Integer> fi = f1;
		while (!fi.isEmpty()) {
			List<Integer> h = new ArrayList<Integer>();
			for (int p : fi) {
				for (int q : S[p]) {
					n[q]--;
					if (n[q] == 0) {
						h.add(q);
					}
				}
			}
			fronts.add(h);
			fi = h;
		}

		int i = 0;
		for (List<Integer> front : fronts) {
			for (int p : front) {
				rank[p] = i;
			}
			i++;
		}

		// To avoid recalculations, save the fronts!
		this.fronts = fronts;

		return fronts;
	}

	protected void calcDistance(List<Integer> list) {
		for (int e : list) {
			dist[e] = 0;
		}

		if (list.size() < 3) {
			return;
		}

		if (m == null) { // initialize the number of objectives
			m = ind[list.get(0)].getObjectives().array().length;
		}

		for (int i = 0; i < m; i++) {
			Collections.sort(list, new DimensionSort(i));
			dist[list.get(0)] = Double.MAX_VALUE;
			dist[list.get(list.size() - 1)] = Double.MAX_VALUE;
			for (int j = 1; j < list.size() - 1; j++) {
				double p = ind[list.get(j - 1)].getObjectives().array()[i];
				double n = ind[list.get(j + 1)].getObjectives().array()[i];
				dist[list.get(j)] += n - p;
			}
		}
	}

	/**
	 * The {@link DistanceSort} sorts {@link Individual}s by the distance in
	 * ascending order.
	 */
	class DistanceSort implements Comparator<Integer> {
		@Override
		public int compare(Integer p, Integer q) {
			Double pv = dist[p];
			Double qv = dist[q];

			return pv.compareTo(qv);
		}
	}

	/**
	 * The {@link DimensionSort} sorts {@link Individual}s by the objective
	 * value {@code d} in ascending order.
	 */
	class DimensionSort implements Comparator<Integer> {

		final int d;

		public DimensionSort(int d) {
			this.d = d;
		}

		@Override
		public int compare(Integer p, Integer q) {
			Double pv = ind[p].getObjectives().array()[d];
			Double qv = ind[q].getObjectives().array()[d];

			return pv.compareTo(qv);
		}
	}

	protected boolean synchronize(Collection<Individual> population) {
		if (population.size() > ind.length) {
			init((int) Math.ceil(population.size() * 1.33));
		}

		Collection<Individual> remove = new LinkedHashSet<Individual>();
		remove.addAll(map.keySet());
		remove.removeAll(population);

		Collection<Individual> add = new LinkedHashSet<Individual>();
		add.addAll(population);
		add.removeAll(map.keySet());

		// unregister remove
		for (Individual e : remove) {
			int i = map.remove(e);
			ind[i] = null;
		}

		// register add
		int i = 0;
		for (Individual e : add) {
			while (ind[i++] != null) {
				;
			}
			ind[--i] = e;
			map.put(e, i);
		}

		return (!remove.isEmpty() || !add.isEmpty());
	}

}
