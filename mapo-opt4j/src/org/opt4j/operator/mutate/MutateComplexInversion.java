package org.opt4j.operator.mutate;

import java.util.Random;

import org.opt4j.genotype.DynamicListGenotype;

import com.google.inject.Inject;

public class MutateComplexInversion implements MutateComplex{

	protected final Random random;
	
	@Inject
	public MutateComplexInversion(Random random) {
		this.random=random;
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public void mutate(DynamicListGenotype<?> genotype, double p) {
		// TODO Auto-generated method stub
		//System.out.println(inversion.toString());
		if(random.nextDouble() < p){
			DynamicListGenotype<Object> inversion = (DynamicListGenotype<Object>) genotype;
			//System.out.println(genotype);
			int alpha=Integer.parseInt(inversion.get(2).toString());
			int idxMutationParent=0; //Mutation index on parent
			int i=4;
			int mutPos=1;
			int sumNuc=0;
			while(alpha!=sumNuc){
				sumNuc+=Integer.parseInt(inversion.get(i).toString());
				i+=3;
				mutPos++;
			}
			idxMutationParent=i+1;
			if(inversion.get(idxMutationParent).equals(")")){
				if(random.nextBoolean()){			
					/*
					 * Select a random position to the left to make the inversion
					 */
					alpha=0;
					int invIdx=0;
					while(invIdx<=0)
							invIdx=random.nextInt(mutPos); //+1
					int j=mutPos;
					i=invIdx;
					while(i>j){
							Object temp= inversion.get((i*3+1));
							if(j==mutPos){
								if(Integer.parseInt(inversion.get((i*3+1)).toString())>= 3 ){
									inversion.set((i*3+1), inversion.get((j*3+1)));
									inversion.set((j*3+1), temp);
								}
							}
							else{
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
							j--;
							i++;
						
					}
					for(i=1; i<mutPos; i++)
						alpha+=Integer.parseInt(inversion.get((i*3+1)).toString());
					inversion.set(2, alpha);
				}
				else{
					/*
					 * Select a random position to the right to make the inversion
					 */
					int invIdx=0;
					while(invIdx<mutPos){
						int pos=(inversion.size()-3)/3;
						invIdx=random.nextInt(pos+1); //+1
					}
					int j=(inversion.size()-3)/3;
					i=mutPos+1;
					while(j>i){
							Object temp= inversion.get((i*3+1));
							//TODO: Change the limit to the real on gene
							if(i==mutPos){
								if(Integer.parseInt(inversion.get((j*3+1)).toString())>=3){
									inversion.set((i*3+1), inversion.get((j*3+1)));
									inversion.set((j*3+1), temp);
								}
							}
							else{
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
							j--;
							i++;
					}
				}
			}
			else if(inversion.get(idxMutationParent).equals(":")){
				/*
				 * Select a random  position to the right to make inversion
				 */
				int invIdx=0;
				while(invIdx<mutPos){
					int pos=(inversion.size()-3)/3;
					invIdx=random.nextInt(pos+1); //+1
				}
				int j=(inversion.size()-3)/3;
				i=mutPos+1;
				while(j>i){
						Object temp= inversion.get((i*3+1));
						//TODO: Change the limit to the real on gene
						if(i==mutPos){
							if(Integer.parseInt(inversion.get((j*3+1)).toString())>=3){
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
						}
						else{
							inversion.set((i*3+1), inversion.get((j*3+1)));
							inversion.set((j*3+1), temp);
						}
						j--;
						i++;
				}
			}
			else if(inversion.get(idxMutationParent).equals(";")){
				/*
				 * Select a random position to the left to make the inversion
				 */
				alpha=0;
				int invIdx=0;
				while(invIdx<=0)
						invIdx=random.nextInt(mutPos); //+1
				int j=mutPos;
				i=invIdx;
				while(i>j){
						Object temp= inversion.get((i*3+1));
						if(j==mutPos){
							if(Integer.parseInt(inversion.get((i*3+1)).toString())>= 3 ){
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
						}
						else{
							inversion.set((i*3+1), inversion.get((j*3+1)));
							inversion.set((j*3+1), temp);
						}
						j--;
						i++;
					
				}
				for(i=1; i<mutPos; i++)
					alpha+=Integer.parseInt(inversion.get((i*3+1)).toString());
				inversion.set(2, alpha);
			}
		}
	}
	
}
