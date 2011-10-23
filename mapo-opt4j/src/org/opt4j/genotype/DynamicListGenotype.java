package org.opt4j.genotype;

import java.util.ArrayList;
import java.util.Collection;

import org.opt4j.core.Genotype;

public class DynamicListGenotype<E> extends ArrayList<E> implements ListGenotype<E> {
	
	public DynamicListGenotype(){
		super();
	}
	
	public DynamicListGenotype(Collection<E> values){
		super(values);
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public <G extends Genotype> G newInstance() {
		try {
			return (G) this.getClass().newInstance();
		} catch (InstantiationException e) {
			throw new RuntimeException(e);
		} catch (IllegalAccessException e) {
			throw new RuntimeException(e);
		}
	}
}
