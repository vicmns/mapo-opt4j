package org.opt4j.operator.crossover;

import java.util.Random;

import org.opt4j.genotype.DynamicListGenotype;
import org.opt4j.optimizer.ea.Pair;

import com.google.inject.Inject;

public class CrossoverComplexOnePoint<G extends DynamicListGenotype<?>> implements Crossover<G>{
	
	protected final Random random;
	
	@Inject
	public CrossoverComplexOnePoint(Random random){
		this.random=random;
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public  Pair<G> crossover(G parent1, G parent2) {
		DynamicListGenotype<Object> o1 = parent1.newInstance();
		DynamicListGenotype<Object> o2 = parent2.newInstance();
		int endIdx=0;
		int idxMparent1=0; //Mutation index on parent 1
		int idxMparent2=0; //Mutation index on parent 2
		int i=4;
		int sumNuc=0;
		//Find the index where mutations start on parents
		while(Integer.parseInt(parent1.get(2).toString())!=sumNuc){
			sumNuc+=Integer.parseInt(parent1.get(i).toString());
			i+=3;
		}
		idxMparent1=i+1;
		i=4;
		sumNuc=0;
		while(Integer.parseInt(parent2.get(2).toString())!=sumNuc){
			sumNuc+=Integer.parseInt(parent2.get(i).toString());
			i+=3;
		}
		idxMparent2=i+1;
		if(random.nextBoolean()){
			//For child 1: Take the left configuration before mutation from parent 1 and right configuration (after) from parent 2
			o1.addAll(parent1.subList(0, idxMparent1+1));
			o1.addAll(parent2.subList(idxMparent2+1, parent2.size()));
			//For child 2: Take the left configuration before mutation from parent 2 and right configuration (after) from parent 1
			o2.addAll(parent2.subList(0, idxMparent2+1));
			o2.addAll(parent1.subList(idxMparent1+1, parent1.size()));
			//Recalculate end Index
			for(i=5; i<o1.size();i+=3){
				endIdx+=Integer.parseInt(o1.get(i-1).toString());
			}
			endIdx+=Integer.parseInt(o1.get(0).toString());
			o1.set(1, endIdx);
			endIdx=0;
			for(i=5;i<o2.size();i+=3){
				endIdx+=Integer.parseInt(o2.get(i-1).toString());
			}
			endIdx+=Integer.parseInt(o2.get(0).toString());
			o2.set(1, endIdx);
		}
		else{
			//For child 1: Take the left configuration before mutation from parent 2 and right configuration (after) from parent 1
			o1.addAll(parent2.subList(0, idxMparent2+1));
			o1.addAll(parent1.subList(idxMparent1+1, parent1.size()));
			//For child 2: Take the left configuration before mutation from parent 1 and right configuration (after) from parent 2
			o2.addAll(parent1.subList(0, idxMparent1+1));
			o2.addAll(parent2.subList(idxMparent2+1, parent2.size()));
			//Recalculate Ending indexes
			for(i=5;i<o1.size();i+=3){
				endIdx+=Integer.parseInt(o1.get(i-1).toString());
			}
			endIdx+=Integer.parseInt(o1.get(0).toString());
			o1.set(1, endIdx);
			endIdx=0;
			for(i=5;i<o2.size();i+=3){
				endIdx+=Integer.parseInt(o2.get(i-1).toString());
			}
			endIdx+=Integer.parseInt(o2.get(0).toString());
			o2.set(1, endIdx);
		}
		Pair<G> offspring = new Pair<G>((G) o1, (G) o2);
		return offspring;
	}
}
