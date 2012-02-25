package org.opt4j.operator.mutate;

import java.util.Random;

import org.mapo.ComplexProblem;
import org.opt4j.genotype.DynamicListGenotype;

import com.google.inject.Inject;

public class MutateComplexInversion implements MutateComplex{

	protected final Random random;
	private int idxMS;
	private int idxME;
	
	@Inject
	public MutateComplexInversion(Random random, ComplexProblem problem) {
		this.random=random;
		this.idxME = problem.getIdxE();
		this.idxMS = problem.getIdxS();
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public void mutate(DynamicListGenotype<?> genotype, double p) {
		if(random.nextDouble() < p){
			DynamicListGenotype<Object> inversion = (DynamicListGenotype<Object>) genotype;
			int alpha=Integer.parseInt(inversion.get(2).toString());
			int idxMutationParent=0; //Mutation index on parent
			int i=4;
			int mutPos=1;
			int sumNuc=0;
			int numInfBulgue=0;
			while(alpha!=sumNuc){
				if(!inversion.get(i-1).equals("{")){
					sumNuc+=Integer.parseInt(inversion.get(i).toString());
				}
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
								if(Integer.parseInt(inversion.get((i*3+1)).toString())>= (this.idxME - this.idxMS + 1) ){
									if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
										numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
										inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
												Integer.parseInt(inversion.get((j*3+1)).toString()) );
									}
									if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
										numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
										inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
												Integer.parseInt(inversion.get((i*3+1)).toString()) );
									}
									inversion.set((i*3+1), inversion.get((j*3+1)));
									inversion.set((j*3+1), temp);
								}
							}
							else{
								if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
									inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
											Integer.parseInt(inversion.get((j*3+1)).toString()) );
								}
								if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
									inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
											Integer.parseInt(inversion.get((i*3+1)).toString()) );
								}
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
							j--;
							i++;
						
					}
					for(i=1; i<mutPos; i++){
						if(!inversion.get(i*3).equals("{")){
							alpha+=Integer.parseInt(inversion.get((i*3+1)).toString());
						}
					}
					int endIdx = 0;
					for(i = mutPos; i < inversion.size() / 3 ; i++){
						if(!inversion.get(i*3).equals("{"))
							endIdx+=Integer.parseInt(inversion.get((i*3+1)).toString());
					}
					inversion.set(2, alpha);
					inversion.set(0, this.idxMS - alpha);
					inversion.set(1, alpha + endIdx + Integer.parseInt(inversion.get(0).toString()) - 1);
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
							Object temp = inversion.get((i*3+1));
							//TODO: Change the limit to the real on gene
							if(i==mutPos){
								if(Integer.parseInt(inversion.get((j*3+1)).toString()) >= (this.idxME - this.idxMS + 1)){
									if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
										numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
										inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
												Integer.parseInt(inversion.get((j*3+1)).toString()) );
									}
									if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
										numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
										inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
												Integer.parseInt(inversion.get((i*3+1)).toString()) );
									}
									inversion.set((i*3+1), inversion.get((j*3+1)));
									inversion.set((j*3+1), temp);
								}
							}
							else{
								if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
									inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
											Integer.parseInt(inversion.get((j*3+1)).toString()) );
								}
								if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
									inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
											Integer.parseInt(inversion.get((i*3+1)).toString()) );
								}
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
							j--;
							i++;
					}
					int endIdx = 0;
					for(i = mutPos; i < inversion.size() / 3 ; i++){
						if(!inversion.get(i*3).equals("{"))
							endIdx+=Integer.parseInt(inversion.get((i*3+1)).toString());
					}
					inversion.set(1, Integer.parseInt(inversion.get(0).toString()) + endIdx + 
							Integer.parseInt(inversion.get(2).toString()) - 1);
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
								if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
									inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
											Integer.parseInt(inversion.get((j*3+1)).toString()) );
								}
								if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
									inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
											Integer.parseInt(inversion.get((i*3+1)).toString()) );
								}
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
						}
						else{
							if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
								numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
								inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
										Integer.parseInt(inversion.get((j*3+1)).toString()) );
							}
							if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
								numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
								inversion.set(1, Integer.parseInt(inversion.get(1).toString()) - numInfBulgue + 
										Integer.parseInt(inversion.get((i*3+1)).toString()) );
							}
							inversion.set((i*3+1), inversion.get((j*3+1)));
							inversion.set((j*3+1), temp);
						}
						j--;
						i++;
				}
				int endIdx = 0;
				for(i = mutPos; i < inversion.size() / 3 ; i++){
					if(!inversion.get(i*3).equals("{"))
						endIdx+=Integer.parseInt(inversion.get((i*3+1)).toString());
				}
				inversion.set(1, Integer.parseInt(inversion.get(0).toString()) + endIdx - 1);
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
				while(i<j){
						Object temp= inversion.get((i*3+1));
						if(j==mutPos){
							if(Integer.parseInt(inversion.get((i*3+1)).toString())>= (this.idxME - this.idxMS + 1) ){
								if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
									inversion.set(0, Integer.parseInt(inversion.get(0).toString()) + numInfBulgue - 
											Integer.parseInt(inversion.get((j*3+1)).toString()) );
								}
								if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
									numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
									inversion.set(0, Integer.parseInt(inversion.get(0).toString()) + numInfBulgue - 
											Integer.parseInt(inversion.get((i*3+1)).toString()) );
								}
								inversion.set((i*3+1), inversion.get((j*3+1)));
								inversion.set((j*3+1), temp);
							}
						}
						else{
							if(inversion.get(i*3).equals("{") && !inversion.get(j*3).equals("{")){
								numInfBulgue = Integer.parseInt(inversion.get((i*3+1)).toString());
								inversion.set(0, Integer.parseInt(inversion.get(0).toString()) + numInfBulgue - 
										Integer.parseInt(inversion.get((j*3+1)).toString()) );
							}
							if(inversion.get(j*3).equals("{") && !inversion.get(i*3).equals("{")){
								numInfBulgue = Integer.parseInt(inversion.get((j*3+1)).toString());
								inversion.set(0, Integer.parseInt(inversion.get(0).toString()) + numInfBulgue - 
										Integer.parseInt(inversion.get((i*3+1)).toString()) );
							}
							inversion.set((i*3+1), inversion.get((j*3+1)));
							inversion.set((j*3+1), temp);
						}
						j--;
						i++;
					
				}
				for(i=1; i<mutPos; i++){
					if(!inversion.get(i*3).equals("{")){
						alpha+=Integer.parseInt(inversion.get((i*3+1)).toString());
					}
				}
				inversion.set(2, alpha);
				inversion.set(0, this.idxME - alpha - Integer.parseInt(inversion.get(mutPos * 3 + 1).toString()) + 1);
			}
		}
	}
	
}
