package org.opt4j.operator.mutate;

import java.util.Random;

import org.opt4j.genotype.DynamicListGenotype;

import com.google.inject.Inject;

public class MutateComplexComplement implements MutateComplex {

	protected Random random;
	
	@Inject
	public MutateComplexComplement(Random random){
		this.random=random;
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public void mutate(DynamicListGenotype<?> genotype, double p) {
		if(random.nextDouble() < p){
			DynamicListGenotype<Object> c= (DynamicListGenotype<Object>) genotype;
			int selectStruct=0;
			int sumNuc=0;
			int i=4;
			int idxMp=0;
			String configuration = "";
			//Find the mutation position
			while(Integer.parseInt(c.get(2).toString())!=sumNuc){
				if(!genotype.get(i-1).equals("{"))
					sumNuc+=Integer.parseInt(c.get(i).toString());
				i+=3;
			}
			idxMp=i+1;
			while(!c.get(selectStruct).equals("[") && !c.get(selectStruct).equals("]"))
				selectStruct=random.nextInt(c.size());
			if(c.get(selectStruct).equals("]")){
				selectStruct-=2;
			}
			configuration = genotype.get(idxMp).toString();
			//Randomly select from which side it will rest
			if(random.nextBoolean() && isConfiguration(c.get(selectStruct-1)) && idxMp!=selectStruct-1){
				//Left side
				int restNuc=random.nextInt(Integer.parseInt(c.get(selectStruct-2).toString()));
				c.set(selectStruct - 2, (Integer.parseInt(c.get(selectStruct-2).toString()) - restNuc));
				c.set(selectStruct + 1, (Integer.parseInt(c.get(selectStruct+1).toString()) + restNuc));
				if(c.get(selectStruct-3).equals("{")){
					if(selectStruct < idxMp){
						c.set(2,Integer.parseInt(c.get(2).toString()) + restNuc );
						if(configuration.equals(";")){
							c.set(0,Integer.parseInt(c.get(0).toString()) - restNuc );
						}
						else
							c.set(1,Integer.parseInt(c.get(1).toString()) + restNuc );
					}
					else
						c.set(1,Integer.parseInt(c.get(1).toString()) + restNuc );
				}
			}
			else if(selectStruct+3<c.size() && isConfiguration(c.get(selectStruct+3)) && idxMp!=selectStruct+5){
				//Right side
				int restNuc=random.nextInt(Integer.parseInt(c.get(selectStruct+4).toString()));
				c.set(selectStruct+4, (Integer.parseInt(c.get(selectStruct+4).toString())-restNuc));
				c.set(selectStruct+1, (Integer.parseInt(c.get(selectStruct+1).toString())+restNuc));
				if(c.get(selectStruct + 3).equals("{")){
					if(selectStruct < idxMp){
						c.set(2,Integer.parseInt(c.get(2).toString()) + restNuc );
						if(configuration.equals(";")){
							c.set(0,Integer.parseInt(c.get(0).toString()) - restNuc );
						}
						else
							c.set(1,Integer.parseInt(c.get(1).toString()) + restNuc );
					}
					else
						c.set(1,Integer.parseInt(c.get(1).toString()) + restNuc );
				}
			}
		}
	}
	
	public boolean isConfiguration(Object item){
		if(item.equals(";"))
			return true;
		if(item.equals(">"))
			return true;
		if(item.equals("}"))
			return true;
		if(item.equals(")"))
			return true;
		else return false;
	}
}
