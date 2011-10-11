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

package org.opt4j.core;

import org.opt4j.core.Individual.State;

/**
 * The {@link IndividualStateListener} receives notifications of an
 * {@link Individual} changing its {@link State}. To add an
 * {@link IndividualStateListener} use the method
 * {@link org.opt4j.start.Opt4JModule#addIndividualStateListener(Class)}.
 * 
 * @author lukasiewycz
 * 
 */
public interface IndividualStateListener {

	/**
	 * Invoked if the {@link Individual} changes its state.
	 * 
	 * @param individual
	 *            the individual that changes the state
	 */
	public void inidividualStateChanged(Individual individual);
}
