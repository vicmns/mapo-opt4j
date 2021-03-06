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

package org.opt4j.common.random;

import java.util.Random;

import org.opt4j.start.Constant;

import com.google.inject.Inject;
import com.google.inject.Singleton;

/**
 * The {@link RandomMersenneTwister} uses an implementation of the mersenne
 * twister random number generator written by {@code David Beaumont}.
 * 
 * @author lukasiewycz
 * 
 */
@SuppressWarnings("serial")
@Singleton
public class RandomMersenneTwister extends MTRandom {

	/**
	 * Constructs a {@link RandomMersenneTwister} with the specified seed.
	 * 
	 * @param seed
	 *            the seed value (using namespace {@link Random})
	 */
	@Inject
	public RandomMersenneTwister(@Constant(value = "seed", namespace = Random.class) long seed) {
		super(seed);
	}

}
