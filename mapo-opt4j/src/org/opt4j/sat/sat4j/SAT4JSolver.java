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

package org.opt4j.sat.sat4j;

import static org.opt4j.sat.Constraint.Operator.EQ;
import static org.opt4j.sat.Constraint.Operator.GE;
import static org.opt4j.sat.Constraint.Operator.LE;

import java.math.BigInteger;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.opt4j.config.annotations.Info;
import org.opt4j.sat.Constraint;
import org.opt4j.sat.Constraint.Operator;
import org.opt4j.sat.ContradictionException;
import org.opt4j.sat.Literal;
import org.opt4j.sat.Model;
import org.opt4j.sat.Order;
import org.opt4j.sat.Solver;
import org.opt4j.sat.TimeoutException;
import org.opt4j.sat.VarOrder;
import org.opt4j.start.Constant;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.minisat.core.LearningStrategy;
import org.sat4j.minisat.core.RestartStrategy;
import org.sat4j.minisat.learning.ClauseOnlyLearning;
import org.sat4j.minisat.learning.FixedLengthLearning;
import org.sat4j.minisat.learning.MiniSATLearning;
import org.sat4j.minisat.orders.PositiveLiteralSelectionStrategy;
import org.sat4j.minisat.orders.VarOrderHeap;
import org.sat4j.minisat.restarts.ArminRestarts;
import org.sat4j.minisat.restarts.LubyRestarts;
import org.sat4j.minisat.restarts.MiniSATRestarts;
import org.sat4j.pb.constraints.CompetResolutionPBMixedHTClauseCardConstrDataStructure;
import org.sat4j.pb.core.PBDataStructureFactory;
import org.sat4j.pb.core.PBSolverResolution;
import org.sat4j.specs.IVec;

import com.google.inject.Inject;

/**
 * The {@link SAT4JSolver} implements a {@link Solver}. It is based on the Java
 * SAT/PB-Solver from <a href="http://www.sat4j.org">Sat4J.org</a>.
 * 
 * @author lukasiewycz
 * 
 */
public class SAT4JSolver implements Solver {

	public enum Learning {
		@Info("Learning only constraints with less or equal number of variables.")
		FIXEDLENGTH, @Info("Use the MiniSAT learning scheme.")
		MINISAT, @Info("Learn only clauses.")
		CLAUSEONLY;
	}

	public enum Restarts {
		@Info("Use MiniSAT restarts.")
		MINISAT, @Info("Use Luby restarts.")
		LUBY, @Info("Use Rapid restarts.")
		RAPID;
	}

	protected final Queue<Constraint> constraints = new ConcurrentLinkedQueue<Constraint>();

	protected final PBSolverResolution solver;

	protected final Map<Object, Integer> variables = new HashMap<Object, Integer>();

	protected int nextVariable = 1;

	/**
	 * Constructs a new {@link SAT4JSolver} with a timeout and the number of
	 * kept learning clauses. Additionally, this constructor allows the
	 * specification of the learning and restart strategy.
	 * 
	 * @param timeout
	 *            timeout in seconds
	 * @param clauseLearningLength
	 *            clauses are learned if they have a smaller or equal number of
	 *            literals per clause
	 * @param learning
	 *            the learning strategy
	 * @param restarts
	 *            the restart strategy
	 */
	@Inject
	public SAT4JSolver(@Constant(value = "timeout", namespace = SAT4JSolver.class) int timeout,
			@Constant(value = "clauseLearningLength", namespace = SAT4JSolver.class) int clauseLearningLength,
			@Constant(value = "learning", namespace = SAT4JSolver.class) Learning learning,
			@Constant(value = "restarts", namespace = SAT4JSolver.class) Restarts restarts) {
		LearningStrategy<PBDataStructureFactory> l = null;
		switch (learning) {
		case FIXEDLENGTH:
			l = new FixedLengthLearning<PBDataStructureFactory>(clauseLearningLength);
			break;
		case MINISAT:
			l = new MiniSATLearning<PBDataStructureFactory>();
			break;
		case CLAUSEONLY:
			l = new ClauseOnlyLearning<PBDataStructureFactory>();
			break;
		}

		RestartStrategy r = null;
		switch (restarts) {
		case MINISAT:
			r = new MiniSATRestarts();
			break;
		case LUBY:
			r = new LubyRestarts();
			break;
		case RAPID:
			r = new ArminRestarts();
			break;
		}
		solver = new PBSolverResolution(l, new CompetResolutionPBMixedHTClauseCardConstrDataStructure(),
				new VarOrderHeap(new PositiveLiteralSelectionStrategy()), r);
		l.setSolver(solver);
		l.setVarActivityListener(solver);
		if (l instanceof MiniSATLearning) {
			((MiniSATLearning<PBDataStructureFactory>) l).setDataStructureFactory(solver.getDSFactory());
		}

		if (timeout <= 0) {
			throw new IllegalArgumentException("Invalid timeout: " + timeout);
		}
		solver.setTimeout(timeout);
		solver.setVerbose(false);

		setNVars(100);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.sat.Solver#addConstraint(org.opt4j.sat.Constraint)
	 */
	@Override
	public void addConstraint(Constraint constraint) {
		constraints.add(constraint);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.opt4j.sat.Solver#solve(org.opt4j.sat.Order)
	 */
	@Override
	public synchronized Model solve(Order order) throws TimeoutException {

		Constraint c = null;
		while ((c = constraints.poll()) != null) {
			addConstraintToSolver(c);
		}

		if (order instanceof VarOrder) {

			VarOrder varOrder = (VarOrder) order;
			VariableOrder o = new VariableOrder();
			solver.setOrder(o);

			for (Entry<Object, Double> entry : varOrder.getActivityEntrySet()) {
				if (variables.containsKey(entry.getKey())) {
					int var = variables.get(entry.getKey());
					o.setVarActivity(var, entry.getValue());
				}
			}

			for (Entry<Object, Boolean> entry : varOrder.getPhaseEntrySet()) {
				if (variables.containsKey(entry.getKey())) {
					int var = variables.get(entry.getKey());
					o.setVarPhase(var, entry.getValue());
				}
			}

			o.setVarInc(varOrder.getVarInc());
			o.setVarDecay(varOrder.getVarDecay());

		}

		try {
			boolean success = solver.isSatisfiable();

			if (success) {
				Model model = new Model();

				for (Entry<Object, Integer> entry : variables.entrySet()) {
					model.set(entry.getKey(), solver.model(entry.getValue()));
				}
				return model;
			} else {
				return null;
			}
		} catch (org.sat4j.specs.TimeoutException e) {
			throw new TimeoutException();
		}
	}

	protected void addConstraintToSolver(Constraint constraint) {

		VecInt lits = toVecInt(constraint.getLiterals());

		IVec<BigInteger> coeffs = new Vec<BigInteger>();

		for (Integer value : constraint.getCoefficients()) {
			coeffs.push(BigInteger.valueOf(value));
		}

		BigInteger d = BigInteger.valueOf(constraint.getRhs());
		Operator operator = constraint.getOperator();

		try {
			if (operator == LE || operator == EQ) {
				solver.addPseudoBoolean(lits, coeffs, false, d);
			}
			if (operator == GE || operator == EQ) {
				solver.addPseudoBoolean(lits, coeffs, true, d);
			}
		} catch (org.sat4j.specs.ContradictionException e) {

			throw new ContradictionException();
		}
	}

	protected VecInt toVecInt(Iterable<Literal> list) {
		VecInt vector = new VecInt();

		for (Literal literal : list) {
			Object var = literal.variable();
			if (!variables.containsKey(var)) {
				variables.put(var, nextVariable++);
				if (variables.size() > solver.nVars()) {
					setNVars(variables.size() * 2);
				}
			}
			boolean phase = literal.phase();

			vector.push(variables.get(var) * (phase ? 1 : -1));
		}

		return vector;
	}

	/**
	 * Sets the instance to n vars.
	 * 
	 * @param n
	 *            the number of vars.
	 */
	protected void setNVars(int n) {
		solver.newVar(n);
	}

}
