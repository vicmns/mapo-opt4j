package org.opt4j.operator.crossover;

import java.util.ArrayList;
import java.util.Random;

import org.mapo.ComplexProblem;
import org.opt4j.genotype.DynamicListGenotype;
import org.opt4j.optimizer.ea.Pair;

import com.google.inject.Inject;

public class CrossoverComplexOnePoint<G extends DynamicListGenotype<?>> implements Crossover<G>{
	
	protected final Random random;
	private int idxMS;
	private int idxME;
	
	@Inject
	public CrossoverComplexOnePoint(Random random, ComplexProblem problem){
		this.random=random;
		this.idxMS=problem.getIdxS();
		this.idxME=problem.getIdxE();
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public  Pair<G> crossover(G parent1, G parent2) {
		DynamicListGenotype<Object> o1 = parent1.newInstance();
		DynamicListGenotype<Object> o2 = parent2.newInstance();
		int idxEC1=0, idxEC2=0, aC1=0, aC2=0, idxSC1=0, idxSC2=0;
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
		//Check which configuration has the parents
		if(parent1.get(idxMparent1).equals(":") && parent2.get(idxMparent2).equals(":")){
			o1.addAll(parent1.subList(idxMparent1-2, idxMparent1+1));
			o1.addAll(parent2.subList(idxMparent2+1, parent2.size()));
			o2.addAll(parent2.subList(idxMparent2-2, idxMparent2+1));
			o2.addAll(parent1.subList(idxMparent1+1, parent1.size()));
			int c1EndIdx=Integer.parseInt(parent2.get(1).toString())-Integer.parseInt(parent2.get(4).toString())+
					Integer.parseInt(parent1.get(4).toString());
			int c2EndIdx=Integer.parseInt(parent1.get(1).toString())-Integer.parseInt(parent1.get(4).toString())+
					Integer.parseInt(parent2.get(4).toString());
			o1.add(0,0);
			o2.add(0,0);
			o1.add(0,c1EndIdx);
			o2.add(0,c2EndIdx);
			o1.add(0,Integer.parseInt(parent2.get(0).toString()));
			o2.add(0,Integer.parseInt(parent1.get(0).toString()));
			System.out.println("OP1");
			if(!o1.get(3).equals(":"))
				System.out.println("I see what wou did there");
			if(!o2.get(3).equals(":"))
				System.out.println("I see what wou did there");
		}
		else if(parent1.get(idxMparent1).equals(":") && parent2.get(idxMparent2).equals(";")){
			int newSIdx=0, newEIdx=0, newAlpha=0;
			o1.addAll(parent1.subList(idxMparent1-2, idxMparent1+1));
			for(int j=idxMparent2-2;j>3;j-=3){
				o1.addAll(parent2.subList(j-3, j));
				newEIdx+=Integer.parseInt(parent2.get(j-2).toString());
			}
			newEIdx+=Integer.parseInt(parent1.get(idxMparent1-1).toString())+Integer.parseInt(parent1.get(0).toString());
			//Repair Genotype
			if(o1.get(0).equals(":") && !o1.get(5).equals("]")){
				o1.set(3, "[");
				o1.set(5, "]");
			}
			repairGenotype(o1);
			o1.add(0,0);
			o1.add(0,newEIdx);
			o1.add(0,parent1.get(0));
			//Child2
			o2.addAll(parent2.subList(idxMparent2-2, idxMparent2+1));
			for(int j=idxMparent1+1;j<parent1.size();j+=3){
				o2.addAll(0,parent1.subList(j, j+3));
				newSIdx+=Integer.parseInt(parent1.get(j+1).toString());
			}
			newAlpha=newSIdx;
			newSIdx=Integer.parseInt(parent2.get(1).toString())-newSIdx-Integer.parseInt(o2.get(o2.size()-2).toString());
			if(newSIdx<0){
				int restToDelete=Math.abs(newSIdx);
				int newNucleo=0;
				int j=1;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o2.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o2.subList(j-1,j+2).clear();
						newAlpha-=Math.abs(lastNucleo);
						newSIdx+=Math.abs(lastNucleo);
					}
					else{
						o2.set(j, newNucleo);
						newAlpha-=Math.abs(restToDelete);
						newSIdx+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			repairGenotype(o2);
			if(o2.get(o2.size()-1).equals(";") &&  !o2.get(o2.size()-4).equals("]")){
				o2.set(o2.size()-6, "[");
				o2.set(o2.size()-4, "]");
			}
			o2.add(0,newAlpha);
			o2.add(0,parent2.get(1));
			o2.add(0,newSIdx);
			System.out.println("OP2");
			if(!o1.get(3).equals(":"))
				System.out.println("I see what wou did there");
			//TODO: OP2 DONE!
		}
		else if(parent1.get(idxMparent1).equals(";") && parent2.get(idxMparent2).equals(":")){
			int newSIdx=0, newEIdx=0, newAlpha=0;
			o1.addAll(parent1.subList(idxMparent1-2, idxMparent1+1));
			for(int j=idxMparent2+1;j<parent2.size();j+=3){
				o1.addAll(0,parent2.subList(j, j+3));
				newSIdx+=Integer.parseInt(parent2.get(j+1).toString());
			}
			newAlpha=newSIdx;
			newSIdx=Integer.parseInt(parent1.get(1).toString())-newSIdx-Integer.parseInt(o1.get(o1.size()-2).toString());
			if(newSIdx<0){
				int restToDelete=Math.abs(newSIdx);
				int newNucleo=0;
				int j=1;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o1.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o1.subList(j-1,j+2).clear();
						newAlpha-=Math.abs(lastNucleo);
						newSIdx+=Math.abs(lastNucleo);
					}
					else{
						o1.set(j, newNucleo);
						newAlpha-=Math.abs(restToDelete);
						newSIdx+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			repairGenotype(o1);
			if(o1.get(o1.size()-1).equals(";") &&  !o1.get(o1.size()-4).equals("]")){
				o1.set(o1.size()-6, "[");
				o1.set(o1.size()-4, "]");
			}
			o1.add(0,newAlpha);
			o1.add(0,parent1.get(1));
			o1.add(0,newSIdx);
			//CHILD 2
			o2.addAll(parent2.subList(idxMparent2-2, idxMparent2+1));
			newEIdx=0;
			for(int j=idxMparent1-2;j>3;j-=3){
				o2.addAll(parent1.subList(j-3, j));
				newEIdx+=Integer.parseInt(parent1.get(j-2).toString());
			}
			newEIdx+=Integer.parseInt(parent2.get(idxMparent2-1).toString())+Integer.parseInt(parent2.get(0).toString());
			if(o2.get(0).equals(":") &&  !o2.get(5).equals("]")){
				o2.set(3, "[");
				o2.set(5, "]");
			}
			repairGenotype(o2);
			o2.add(0,0);
			o2.add(0,newEIdx);
			o2.add(0,parent2.get(0));
			System.out.println("OP3");
			if(!o2.get(3).equals(":"))
				System.out.println("I see what wou did there");
			//TODO: OP3 DONE!
		}
		else if(parent1.get(idxMparent1).equals(";") && parent2.get(idxMparent2).equals(";")){
			int newSIdx=0, newEIdx=0, newAlpha=0;
			o1.addAll(parent2.subList(0, idxMparent2-2));
			o2.addAll(parent1.subList(0, idxMparent1-2));
			o1.addAll(parent1.subList(idxMparent1-2, idxMparent1+1));
			o2.addAll(parent2.subList(idxMparent2-2, idxMparent2+1));
			//Check if the genotype of child 1 is correct
			newAlpha=Integer.parseInt(o1.get(2).toString());
			newEIdx=newAlpha+Integer.parseInt(o1.get(o1.size()-2).toString());
			newSIdx=newAlpha;
			//Check if the genotype of child 2 is correct
			newAlpha=Integer.parseInt(o2.get(2).toString());
			newEIdx=newAlpha+Integer.parseInt(o2.get(o2.size()-2).toString());
			newSIdx=newAlpha;
			System.out.println("OP4");
			//TODO: OP4 DONE!
		}
		else if(parent1.get(idxMparent1).equals(":") && parent2.get(idxMparent2).equals(")")){
			//Left Overhang
			o1.addAll(parent1.subList(0, idxMparent1+1));
			o1.addAll(parent2.subList(idxMparent2+1, parent2.size()));
			o2.addAll(parent2.subList(0, idxMparent2+1));
			o2.addAll(parent1.subList(idxMparent1+1, parent1.size()));
			//Recalculate end indexes
			for(i=5; i<o1.size();i+=3){
				idxEC1+=Integer.parseInt(o1.get(i-1).toString());
			}
			idxEC1+=Integer.parseInt(o1.get(0).toString());
			o1.set(1, idxEC1);
			for(i=5;i<o2.size();i+=3){
				idxEC2+=Integer.parseInt(o2.get(i-1).toString());
			}
			idxEC2+=Integer.parseInt(o2.get(0).toString());
			o2.set(1, idxEC2);
			System.out.println("OP5");
			if(!o1.get(3).equals(":"))
				System.out.println("I see what wou did there");
			//OP5 DONE!
		}
		else if(parent1.get(idxMparent1).equals(")") && parent2.get(idxMparent2).equals(":")){
			//Left Overhang
			o1.addAll(parent1.subList(0, idxMparent1+1));
			o1.addAll(parent2.subList(idxMparent2+1, parent2.size()));
			o2.addAll(parent2.subList(0, idxMparent2+1));
			o2.addAll(parent1.subList(idxMparent1+1, parent1.size()));
			//Recalculate end indexes
			for(i=5; i<o1.size();i+=3){
				idxEC1+=Integer.parseInt(o1.get(i-1).toString());
			}
			idxEC1+=Integer.parseInt(o1.get(0).toString());
			o1.set(1, idxEC1);
			for(i=5;i<o2.size();i+=3){
				idxEC2+=Integer.parseInt(o2.get(i-1).toString());
			}
			idxEC2+=Integer.parseInt(o2.get(0).toString());
			o2.set(1, idxEC2);
			System.out.println("OP6");
			if(!o2.get(3).equals(":"))
				System.out.println("I see what wou did there");
			//TODO: OP6 DONE!
		}
		else if(parent1.get(idxMparent1).equals(";") && parent2.get(idxMparent2).equals(")")){
			o1.addAll(parent2.subList(0, idxMparent2-2));
			o1.addAll(parent1.subList(idxMparent1-2, idxMparent1+1));
			o2.addAll(parent1.subList(0, idxMparent1-2));
			o2.addAll(parent2.subList(idxMparent2-2,parent2.size()));
			//...,===>(#)...[#]<===
			//Recalculate Start Indexes
			/*for(i=4; i<idxMparent2-1;i+=3){
				aC1+=Integer.parseInt(parent2.get(i).toString());
			}*/
			aC1=Integer.parseInt(o1.get(2).toString());
			idxEC1+=aC1+Integer.parseInt(o1.get(o1.size()-2).toString());
			/*for(i=4; i<idxMparent1-1;i+=3){
				aC2+=Integer.parseInt(parent1.get(i).toString());
			}*/
			aC2=Integer.parseInt(o2.get(2).toString());
			//END index for child 2
			for(i=idxMparent2-1; i<parent2.size();i+=3){
				idxEC2+=Integer.parseInt(parent2.get(i).toString());
			}
			idxEC2+=aC2;
			//Check if the genotype of child 1 is correct
			idxSC1=idxME-aC1-Integer.parseInt(o1.get(o1.size()-2).toString());
			if(idxSC1<0){
				int restToDelete=Math.abs(idxSC1);
				int newNucleo=0;
				int j=4;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o1.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o1.subList(j-1,j+2).clear();
						aC1-=Math.abs(lastNucleo);
						idxSC1+=Math.abs(lastNucleo);
					}
					else{
						o1.set(j, newNucleo);
						aC1-=Math.abs(restToDelete);
						idxSC1+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			repairGenotypeWithIndexes(o1);
			//Check if the genotype of child 2 is correct
			idxSC2=this.idxMS-aC2;
			if(idxSC2<0){
				int restToDelete=Math.abs(idxSC2);
				int newNucleo=0;
				int j=4;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o2.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o2.subList(j-1,j+2).clear();
						aC2-=Math.abs(lastNucleo);
						idxSC2+=Math.abs(lastNucleo);
					}
					else{
						o2.set(j, newNucleo);
						aC2-=Math.abs(restToDelete);
						idxSC2+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			repairGenotypeWithIndexes(o2);
			//Add indexes
			o1.set(2,aC1);
			o1.set(1,this.idxME);
			o1.set(0,idxSC1);
			o2.set(2,aC2);
			o2.set(1,idxEC2+idxSC2);
			o2.set(0,idxSC2);		
			System.out.println("OP7");
			//TODO: OP7 DONE!
		}
		else if(parent1.get(idxMparent1).equals(")") && parent2.get(idxMparent2).equals(";")){
			o1.addAll(parent2.subList(0, idxMparent2-2));
			o1.addAll(parent1.subList(idxMparent1-2, parent1.size()));
			o2.addAll(parent1.subList(0, idxMparent1-2));
			o2.addAll(parent2.subList(idxMparent2-2,idxMparent2+1));
			//Recalculate Start Indexes
			/*for(i=4; i<idxMparent2-1;i+=3){
				aC1+=Integer.parseInt(parent2.get(i).toString());
			}*/
			aC1=Integer.parseInt(o1.get(2).toString());
			/*for(i=4; i<idxMparent1-1;i+=3){
				aC2+=Integer.parseInt(parent1.get(i).toString());
			}*/
			aC2=Integer.parseInt(o2.get(2).toString());
			//END index for child 1
			for(i=idxMparent1-1; i<parent1.size();i+=3){
				idxEC1+=Integer.parseInt(parent1.get(i).toString());
			}
			idxEC1+=aC1+Integer.parseInt(parent1.get(idxMparent1-1).toString());
			idxEC2+=aC2+Integer.parseInt(o2.get(o2.size()-2).toString());
			//Check if the genotype of child 1 is correct
			idxSC1=this.idxMS-aC1;
			if(idxSC1<0){
				int restToDelete=Math.abs(idxSC1);
				int newNucleo=0;
				int j=4;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o1.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o1.subList(j-1,j+2).clear();
						aC1-=Math.abs(lastNucleo);
						idxSC1+=Math.abs(lastNucleo);
					}
					else{
						o1.set(j, newNucleo);
						aC1-=Math.abs(restToDelete);
						idxSC1+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			repairGenotypeWithIndexes(o1);
			//Check if the genotype of child 2 is correct
			idxSC2=idxME-aC2-Integer.parseInt(o2.get(o2.size()-2).toString());;
			if(idxSC2<0){
				int restToDelete=Math.abs(idxSC2);
				int newNucleo=0;
				int j=4;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o2.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o2.subList(j-1,j+2).clear();
						aC2-=Math.abs(lastNucleo);
						idxSC2+=Math.abs(lastNucleo);
					}
					else{
						o2.set(j, newNucleo);
						aC2-=Math.abs(restToDelete);
						idxSC2+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			repairGenotypeWithIndexes(o2);
			//Add indexes
			o1.set(2,aC1);
			o1.set(1,idxEC1+idxSC1);
			o1.set(0,idxSC1);
			o2.set(2,aC2);
			o2.set(1,this.idxME);
			o2.set(0,idxSC2);
			System.out.println("OP8");
			//TODO: OP8 DONE! (NOPE CHUCK TESTA!)
		}		
		else if(parent1.get(idxMparent1).equals(")") && parent2.get(idxMparent2).equals(")")){
			if(random.nextBoolean()){
				/*For child 1: Take the left configuration from start to Mutation
				 * From parent 1 and right configuration (after) from parent 2
				 */
				o1.addAll(parent1.subList(0, idxMparent1+1));
				o1.addAll(parent2.subList(idxMparent2+1, parent2.size()));
				//For child 2: Take the left configuration before mutation from parent 2 and right configuration (after) from parent 1
				o2.addAll(parent2.subList(0, idxMparent2+1));
				o2.addAll(parent1.subList(idxMparent1+1, parent1.size()));
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
				o1.addAll(parent2.subList(0, idxMparent2+1));
				o1.addAll(parent1.subList(idxMparent1+1, parent1.size()));
				//For child 2: Take the left configuration before mutation from parent 1 and right configuration (after) from parent 2
				o2.addAll(parent1.subList(0, idxMparent1+1));
				o2.addAll(parent2.subList(idxMparent2+1, parent2.size()));
				//Recalculate Ending indexes
				int endIdx=0;
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
		}
		else{
			System.out.println("Algo indesperado pasó aquí");
		}
		Pair<G> offspring = new Pair<G>((G) o1, (G) o2);
		return offspring;
	}
	
	private void repairGenotype(ArrayList<Object> c){
		if(c.get(0).equals(";")){
			c.set(0,":");
			c.set(2, ":");
		}
		if(c.get(0).equals(";")){
			c.set(0,":");
			c.set(2, ":");
		}
		if(c.get(0).equals(":") && !c.get(3).equals("[")){
			c.set(0,"[");
			c.set(2, "]");
		}
		if(c.get(3).equals("[") && c.get(0).equals("[")){
			int tempNuc=Integer.parseInt(c.get(1).toString());
			c.set(4, (Integer.parseInt(c.get(4).toString())+tempNuc));
			c.subList(0,3).clear();
		}
		if(c.get(c.size()-1).equals(":") ){
			c.set(c.size()-3,";");
			c.set(c.size()-1, ";");
		}
		if(!c.get(c.size()-4).equals("]") && c.get(c.size()-1).equals(";")){
			c.set(c.size()-3,"[");
			c.set(c.size()-1, "]");
		}
		if(c.get(c.size()-4).equals("]") && c.get(c.size()-1).equals("]")){
			int tempNuc=Integer.parseInt(c.get(c.size()-2).toString());
			c.set(c.size()-5, (Integer.parseInt(c.get(c.size()-5).toString())+tempNuc));
			c.subList(c.size()-3, c.size()).clear();
			
		}
	}
	
	private void repairGenotypeWithIndexes(ArrayList<Object> c){
		if(c.get(3).equals(";")){
			c.set(3,":");
			c.set(5, ":");
		}
		if(c.get(3).equals(";")){
			c.set(3,":");
			c.set(5, ":");
		}
		if(!c.get(6).equals("[") && c.get(3).equals(":")){
			c.set(3,"[");
			c.set(5, "]");
		}
		if(c.get(6).equals("[") && c.get(3).equals("[")){
			int tempNuc=Integer.parseInt(c.get(4).toString());
			c.set(7, (Integer.parseInt(c.get(7).toString())+tempNuc));
			c.subList(3,6).clear();
		}
		if(c.get(c.size()-1).equals(":") ){
			c.set(c.size()-3,";");
			c.set(c.size()-1, ";");
		}
		if(!c.get(c.size()-4).equals("]") && c.get(c.size()-1).equals(";")){
			c.set(c.size()-3,"[");
			c.set(c.size()-1, "]");
		}
		if(c.get(c.size()-4).equals("]") && c.get(c.size()-1).equals("]")){
			int tempNuc=Integer.parseInt(c.get(c.size()-2).toString());
			c.set(c.size()-5, (Integer.parseInt(c.get(c.size()-5).toString())+tempNuc));
			c.subList(c.size()-3, c.size()).clear();
			
		}
	}
}
