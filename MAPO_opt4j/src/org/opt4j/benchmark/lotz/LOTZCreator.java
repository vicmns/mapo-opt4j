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

package org.opt4j.benchmark.lotz;

import java.util.Random;

import org.opt4j.benchmark.BinaryString;
import org.opt4j.common.random.Rand;
import org.opt4j.core.problem.Creator;
import org.opt4j.start.Constant;

import com.google.inject.Inject;

/**
 * The {@link LOTZCreator}.
 * 
 * @author lukasiewycz
 * 
 */
public class LOTZCreator implements Creator<BinaryString> {

	protected final int size;

	protected final Random random;

	/**
	 * Constructs a {@link LOTZCreator}.
	 * 
	 * @param random
	 *            the random number generator
	 * @param size
	 *            the size of the string
	 */
	@Inject
	public LOTZCreator(Rand random, @Constant(value = "size", namespace = LOTZCreator.class) int size) {
		super();
		this.random = random;
		this.size = size;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.core.problem.Creator#create()
	 */
	@Override
	public BinaryString create() {
		BinaryString string = new BinaryString();
		string.init(random, size);

		return string;
	}
}
