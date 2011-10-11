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

package org.opt4j.operator.algebra;

/**
 * The {@link Pow} performs an exponentiation of two {@link Term}s.
 * 
 * @author lukasiewycz
 * 
 */
public class Pow implements Term {

	protected final Term base;

	protected final Term exponent;

	/**
	 * Constructs a {@link Pow} term.
	 * 
	 * @param base
	 *            the base
	 * @param exponent
	 *            the exponent
	 */
	public Pow(final Term base, final Term exponent) {
		super();
		this.base = base;
		this.exponent = exponent;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.operator.algebra.Term#calculate(double[])
	 */
	@Override
	public double calculate(double... values) {
		return Math.pow(base.calculate(values), exponent.calculate(values));
	}

}
