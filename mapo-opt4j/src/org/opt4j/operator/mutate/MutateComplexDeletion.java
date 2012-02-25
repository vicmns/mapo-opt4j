package org.opt4j.operator.mutate;

import java.util.Random;

import org.opt4j.genotype.DynamicListGenotype;

import com.google.inject.Inject;

public class MutateComplexDeletion implements MutateComplex{

	Random random;
	
	@Inject
	MutateComplexDeletion(Random random){
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
			int numInfBulgue=0;
			//Find the mutation position
			while(Integer.parseInt(c.get(2).toString())!=sumNuc){
				if(!genotype.get(i-1).equals("{"))
					sumNuc+=Integer.parseInt(c.get(i).toString());
				i+=3;
			}
			idxMp=i+1;
			/**
			 * Verify if complex configuration is greater than 3 complexes i.e.
			 * i,j,a,:,3,:,[,4,],;,5,;
			 */
			if((c.size())/3>5){
				int a=0, b=0;
				boolean side=false, strucDeleted=false;
				do{
					selectStruct=random.nextInt((c.size()/3));
					selectStruct=selectStruct*3+2;
				}while(!isConfiguration(c.get(selectStruct)) || idxMp==selectStruct);
				/**
				 * Select a random position (left or right) to delete
				 */
				if(random.nextBoolean() && selectStruct+1<c.size()){
					//Right
					a=Integer.parseInt(c.get(selectStruct-1).toString());
					b=Integer.parseInt(c.get(selectStruct+2).toString());
					if(c.get(selectStruct).equals("}")){
						numInfBulgue = a;
					}
					side=true;				
					strucDeleted=true;
					
				}
				else if(selectStruct-3!=idxMp && selectStruct-5>2 && selectStruct+1<c.size()){
					//Left
					a=Integer.parseInt(c.get(selectStruct-4).toString());
					b=Integer.parseInt(c.get(selectStruct-1).toString());
					if(c.get(selectStruct).equals("}")){
						numInfBulgue = b;
					}
					strucDeleted=true;
				}
				if(strucDeleted){
					//Recalculate indexes
					int endIdx = Integer.parseInt(c.get(1).toString());
					int alpha = Integer.parseInt(c.get(2).toString());
					int startIdx = Integer.parseInt(c.get(0).toString());
					if(selectStruct>idxMp){
						//Recalculate end index
						endIdx-=(a+b-numInfBulgue);
						c.set(1, endIdx);
					}
					else{
						//Recalculate start index and alpha index
						alpha-=(a+b-numInfBulgue);
						c.set(2,alpha);
						startIdx+=(a+b-numInfBulgue);
						c.set(0, startIdx);
						
					}
					if(side)
						c.subList(selectStruct-2, selectStruct+4).clear();
					else
						c.subList(selectStruct-5, selectStruct+1).clear();
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
