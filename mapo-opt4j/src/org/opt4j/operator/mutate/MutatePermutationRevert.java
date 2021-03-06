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

package org.opt4j.operator.mutate;

import java.util.Collections;
import java.util.Random;

import org.opt4j.common.random.Rand;
import org.opt4j.genotype.PermutationGenotype;

import com.google.inject.Inject;

/**
 * <p>
 * Mutate for the {@link PermutationGenotype}. With a given mutation rate two
 * elements are selected and the list between is inverted.
 * </p>
 * 
 * <p>
 * Given a permutation {@code 1 2 3 4 5 6 7 8}, this might result in {@code 1 2
 * 7 6 5 3 4}.
 * </p>
 * 
 * @author lukasiewycz
 * 
 */
public class MutatePermutationRevert implements MutatePermutation {

	protected final Random random;

	/**
	 * Constructs a new {@link MutatePermutation} with the given mutation rate.
	 * 
	 * @param random
	 *            the random number generator
	 */
	@Inject
	public MutatePermutationRevert(Rand random) {
		this.random = random;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.opt4j.operator.mutate.Mutate#mutate(org.opt4j.core.problem.Genotype,
	 * double)
	 */
	@Override
	public void mutate(PermutationGenotype<?> genotype, double p) {
		int size = genotype.size();

		if (size > 1) {
			for (int a = 0; a < size - 1; a++) {
				if (random.nextDouble() < p) {
					int b;
					do {
						b = a + random.nextInt(size - a);
					} while (b == a);

					while (a < b) {
						Collections.swap(genotype, a, b);
						a++;
						b--;
					}
				}
			}
		}
	}

}