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

package org.opt4j.benchmark.zdt;

import static java.lang.Math.pow;

/**
 * Function ZDT 2.
 * 
 * @author lukasiewycz
 * 
 */
public class ZDT2 extends ZDT1 {

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.benchmark.zdt.ZDT1#h(double, double)
	 */
	@Override
	protected double h(double f1, double g) {
		double h = 1.0 - pow(f1 / g, 2.0);
		return h;
	}

}
