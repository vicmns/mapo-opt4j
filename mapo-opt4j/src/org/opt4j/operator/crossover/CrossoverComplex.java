package org.opt4j.operator.crossover;

import java.util.Random;

import org.opt4j.common.random.Rand;
import org.opt4j.genotype.DynamicListGenotype;
import org.opt4j.optimizer.ea.Pair;

public class CrossoverComplex<G extends DynamicListGenotype<?>> implements Crossover<G> {
	protected final Random random;
	
	public CrossoverComplex(Rand random) {
		this.random=random;
	}
	
	public  Pair<G> crossover(G p1, G p2) {
		DynamicListGenotype<Object> o1 = p1.newInstance();
		DynamicListGenotype<Object> o2 = p2.newInstance();
		
		int idxMp1=0; //Mutation index on parent 1
		int idxMp2=0; //Mutation index on parent 2
		int i=4;
		int sumNuc=0;
		//Find the index where mutations start on parents
		while(Integer.parseInt(p1.get(2).toString())!=sumNuc){
			sumNuc+=Integer.parseInt(p1.get(i).toString());
			i+=3;
		}
		idxMp1=i+1;
		i=4;
		sumNuc=0;
		while(Integer.parseInt(p2.get(2).toString())!=sumNuc){
			sumNuc+=Integer.parseInt(p2.get(i).toString());
			i+=3;
		}
		idxMp2=i+1;
		if(random.nextBoolean()){
			//For child 1: Take the left configuration before mutation from parent 1 and right configuration (after) from parent 2
			o1.addAll(p1.subList(0, idxMp1+1));
			o1.addAll(p2.subList(idxMp2+1, p2.size()));
			//For child 2: Take the left configuration before mutation from parent 2 and right configuration (after) from parent 1
			o2.addAll(p2.subList(0, idxMp2+1));
			o2.addAll(p1.subList(idxMp1+1, p1.size()));
			//Recalculate end Index
			int endIdx=0;
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
			o1.addAll(p2.subList(0, idxMp2+1));
			o1.addAll(p1.subList(idxMp1+1, p1.size()));
			//For child 2: Take the left configuration before mutation from parent 1 and right configuration (after) from parent 2
			o2.addAll(p1.subList(0, idxMp1+1));
			o2.addAll(p2.subList(idxMp2+1, p2.size()));
			//Recalculate Ending indexes
			int endIdx=0;
			for(i=5;i<o1.size();i+=3){
				endIdx+=Integer.parseInt(o1.get(0).toString());
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
		@SuppressWarnings("unchecked")
		Pair<G> offspring = new Pair<G>((G) o1, (G) o2);
		return offspring;
	}
}
