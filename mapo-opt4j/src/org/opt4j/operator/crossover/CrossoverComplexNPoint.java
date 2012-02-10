package org.opt4j.operator.crossover;

import java.util.ArrayList;
import java.util.Random;

import org.mapo.ComplexProblem;
import org.opt4j.genotype.DynamicListGenotype;
import org.opt4j.optimizer.ea.Pair;

public class CrossoverComplexNPoint<G extends DynamicListGenotype<?>> implements Crossover<G> {

	protected final Random random;
	private int idxMS;
	private int idxME;
	
	public CrossoverComplexNPoint(Random random, ComplexProblem problem){
		this.random=random;
		this.idxMS=problem.getIdxS();
		this.idxME=problem.getIdxE();
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public  Pair<G> crossover(G parent1, G parent2) {
		DynamicListGenotype<Object> o1 = parent1.newInstance();
		DynamicListGenotype<Object> o2 = parent2.newInstance();
		
		int aC1=0;
		int aC2=0;
		int idxMp1=0; //Mutation index on parent 1
		int idxMp2=0; //Mutation index on parent 2
		int idxMc1=0; //Mutation index on child 1
		int idxMc2=0; //Mutation index on child 2
		int idxEC1=0, idxEC2=0;
		int idxSC1=0, idxSC2=0;
		int i=4; //Initial Index on Individual where the complex configuration start
		/*
		 * [i,j,a,;,#,;[,#,],...]
		 * [0,1,2,3,4,...]
		 */
		int sumNuc=0;
		//Find the index where mutations start on parents
		while(Integer.parseInt(parent1.get(2).toString())!=sumNuc){
			sumNuc+=Integer.parseInt(parent1.get(i).toString());
			i+=3;
		}
		idxMp1=i+1;
		/*
		 * [i,j,a,:,#,:[,#,],(,#,)<==,...]
		 * [i,j,a,*:,#,:[,#,],(,#,),...]
		 * [i,j,a,:#,:,[,#,],(,#,),...,[,#,],;,#,;<==]
		 */
		int idxMpT1=idxMp1+3;
		/*
		 * [i,j,a,;,#,;[,#,],(,#,),[,#,]<==...]
		 * [i,j,a,:,#,:<==[,#,],(,#,),...]
		 * [i,j,a,:#,:[,#,],(,#,),...,[,#,],;,#,;*]
		 */
		i=4;
		sumNuc=0;
		while(Integer.parseInt(parent2.get(2).toString())!=sumNuc){
			sumNuc+=Integer.parseInt(parent2.get(i).toString());
			i+=3;
		}
		idxMp2=i+1;
		int idxMpT2=idxMp2+3;
		/*
		 * Crossover on different cases of representation
		 */
		if(parent1.get(idxMp1).equals(":") && parent2.get(idxMp2).equals(":")){
			//System.out.println("OP1");
			if(random.nextBoolean()){
				o1.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
				o2.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
				
			}
			else{
				o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
				o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
			}
			//Crossover on the right side
			while(idxMpT1 <= parent1.size()-4 && idxMpT2 <= parent2.size()-4){
				if(random.nextBoolean()){
					o1.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					idxEC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					idxEC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					idxEC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					idxEC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1+=3;
				idxMpT2+=3;
			}
			//Align the terminal, then crossover
			if(random.nextBoolean()){
				o1.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				o2.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				idxEC1+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
				idxEC2+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
			}
			else{
				o1.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				o2.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				idxEC1+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
				idxEC2+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
			}
			if((idxMpT2-parent2.size())<(idxMpT1-parent1.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o1.addAll(o1.size()-3,parent2.subList(idxMpT2-2, parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC1+=Integer.parseInt(parent2.get(j).toString());					
				}
				else{
					o2.addAll(o2.size()-3,parent2.subList(idxMpT2-2,parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC2+=Integer.parseInt(parent2.get(j).toString());		
				}
			}
			if((idxMpT1-parent1.size())<(idxMpT2-parent2.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o2.addAll(o2.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC2+=Integer.parseInt(parent1.get(j).toString());	
				}
				else{
					o1.addAll(o1.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC1+=Integer.parseInt(parent1.get(j).toString());	
				}
			}
			idxEC1+=Integer.parseInt(o1.get(1).toString()) - 1;
			idxEC2+=Integer.parseInt(o2.get(1).toString()) - 1;
			//Add indexes
			o1.add(0,0);
			o1.add(0,idxEC1+this.idxMS);
			o1.add(0,this.idxMS);
			o2.add(0,0);
			o2.add(0,idxEC2+this.idxMS);
			o2.add(0,this.idxMS);
			//Repair the genotype
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP1 DONE!
		}
		else if(parent1.get(idxMp1).equals(":") && parent2.get(idxMp2).equals(";")){
			/*
			 * [i,j,a,:,#,:*[,#,],(,#,),...]
			 * [i,j,a,:#,:[,#,],(,#,),...,;,#,;*]
			 */
			//System.out.println("OP2");
			o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
			o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
			idxMp2-=3;
			idxMpT2-=6;
			//Crossover on both sides
			while(idxMpT1 < parent1.size()-1 && idxMpT2 > 5){
				if(random.nextBoolean()){
					o1.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					idxEC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					aC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					idxEC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					aC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1+=3;
				idxMpT2-=3;
			}
			if(random.nextBoolean()){
				//CrossOver terminals blocks
				for(int j=idxMpT1-2; j<parent1.size(); j+=3){
					o2.addAll(0,parent1.subList(j, j+3));
					aC2+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j>=3; j-=3){
					o1.addAll(parent2.subList(j, j+3));
					idxEC1+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			else{
				for(int j=idxMpT2-2; j>=3; j-=3){
					o2.addAll(0,parent2.subList(j, j+3));
					aC2+=Integer.parseInt(parent2.get(j+1).toString());
				}
				for(int j=idxMpT1-2; j<parent1.size(); j+=3){
					o1.addAll(parent1.subList(j, j+3));
					idxEC1+=Integer.parseInt(parent1.get(j+1).toString());
				}
			}
			//idxEC2+=Integer.parseInt(o2.get(o2.size()-2).toString());
			//Check if the genotype of child 2 is correct
			idxSC2=this.idxME-aC2-Integer.parseInt(o2.get(o2.size()-2).toString()) + 1;
			if(idxSC2<0){
				int restToDelete=Math.abs(idxSC2);
				int newNucleo=0;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o2.get(1).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o2.subList(0,3).clear();
						aC2-=Math.abs(lastNucleo);
						idxSC2+=Math.abs(lastNucleo);
					}
					else{
						o2.set(1, newNucleo);
						aC2-=Math.abs(restToDelete);
						idxSC2+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			idxEC1+=Integer.parseInt(o1.get(1).toString()) -1;
			//aC2=this.idxME-Integer.parseInt(o2.get(o2.size()-2).toString());
			o1.add(0,0);
			o1.add(0,idxEC1+this.idxMS);
			o1.add(0,this.idxMS);
			o2.add(0,aC2);
			o2.add(0,this.idxME);
			o2.add(0,idxSC2);
			//Repair Genotypes
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP2 DONE!
		}
		else if(parent1.get(idxMp1).equals(";") && parent2.get(idxMp2).equals(":")){
			/*
			 * [i,j,a,:,#,:*[,#,],(,#,),...]
			 * [i,j,a,:#,:[,#,],(,#,),...,;,#,;*]
			 */
			//System.out.println("OP3");
			o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
			o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
			idxMp1-=3;
			idxMpT1-=6;
			//Crossover on both sides
			while(idxMpT1 > 5 && idxMpT2 < parent2.size()-1){
				if(random.nextBoolean()){
					o1.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					aC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					idxEC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					aC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					idxEC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1-=3;
				idxMpT2+=3;
			}
			//Cross the rest
			if(random.nextBoolean()){
				//Cross end blocks
				for(int j=idxMpT1-2; j>=3; j-=3){
					o2.addAll(parent1.subList(j, j+3));
					idxEC2+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j<parent2.size(); j+=3){
					o1.addAll(0,parent2.subList(j, j+3));
					aC1+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			else{
				//Stay the same way
				for(int j=idxMpT1-2; j>=3; j-=3){
					o1.addAll(0,parent1.subList(j, j+3));
					aC1+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j<parent2.size(); j+=3){
					o2.addAll(parent2.subList(j, j+3));
					idxEC2+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			//idxEC1+=Integer.parseInt(o1.get(o1.size()-2).toString());
			//Check if the genotype of child 1 is correct
			idxSC1=this.idxME-aC1-Integer.parseInt(o1.get(o1.size()-2).toString()) + 1;
			if(idxSC1<0){
				int restToDelete=Math.abs(idxSC1);
				int newNucleo=0;
				int j=1;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o1.get(1).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						o1.subList(0,3).clear();
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
			idxEC2+=Integer.parseInt(o2.get(1).toString()) - 1;
			o1.add(0,aC1);
			o1.add(0,this.idxME);
			o1.add(0,idxSC1);
			o2.add(0,0);
			o2.add(0,idxEC2+this.idxMS);
			o2.add(0,this.idxMS);
			//Repair Genotypes
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP3 Done!
		}
		else if(parent1.get(idxMp1).equals(";") && parent2.get(idxMp2).equals(";")){
			//System.out.println("OP4");
			idxMpT1-=3;
			idxMpT2-=3;
			if(random.nextBoolean()){
				o1.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
				o2.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
				idxMpT1-=3;
				idxMpT2-=3;
				
			}
			else{
				o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
				o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
				idxMpT1-=3;
				idxMpT2-=3;
			}
			//Next crossover on the left side of complexes
			while(idxMpT1>5 && idxMpT2>5){
				if(random.nextBoolean()){
					o1.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					aC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					aC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					aC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					aC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1-=3;
				idxMpT2-=3;
			}
			//CrossOver the rest of configurations
			if(random.nextBoolean()){
				//Cross initiators blocks
				for(int j=idxMpT1-2; j>=3; j-=3){
					o2.addAll(0,parent1.subList(j, j+3));
					aC2+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j>=3; j-=3){
					o1.addAll(0,parent2.subList(j,j+3));
					aC1+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			else{
				for(int j=idxMpT1-2; j>=3; j-=3){
					o1.addAll(0,parent1.subList(j, j+3));
					aC1+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j>=3; j-=3){
					o2.addAll(0,parent2.subList(j, j+3));
					aC2+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			//idxEC1+=aC1+Integer.parseInt(o1.get(o1.size()-2).toString());
			//idxEC2+=aC2+Integer.parseInt(o2.get(o2.size()-2).toString());
			//Check if the genotype of child 1 is correct
			idxSC1=this.idxME-aC1-Integer.parseInt(o1.get(o1.size()-2).toString()) + 1;
			if(idxSC1<0){
				int restToDelete=Math.abs(idxSC1);
				int newNucleo=0;
				int j=1;
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
			//Check if the genotype of child 2 is correct
			idxSC2=this.idxME-aC2-Integer.parseInt(o2.get(o2.size()-2).toString()) + 1;
			if(idxSC2<0){
				int restToDelete=Math.abs(idxSC2);
				int newNucleo=0;
				int j=1;
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
			//Add indexes
			o1.add(0,aC1);
			o1.add(0,this.idxME);
			o1.add(0,idxSC1);
			o2.add(0,aC2);
			o2.add(0,this.idxME);
			o2.add(0,idxSC2);
			//Repair the genotypes
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP4 DONE!
		}
		else if(parent1.get(idxMp1).equals(":") && parent2.get(idxMp2).equals(")")){
			//System.out.println("OP5");
			o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
			o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
			//Crossover on the right side
			while(idxMpT1 <= parent1.size()-4 && idxMpT2 <= parent2.size()-4){
				if(random.nextBoolean()){
					o1.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					idxEC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					idxEC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					idxEC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					idxEC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1+=3;
				idxMpT2+=3;
			}
			
			//Align the terminal, then crossover
			if(random.nextBoolean()){
				o1.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				o2.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				idxEC1+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
				idxEC2+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
			}
			else{
				o1.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				o2.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				idxEC1+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
				idxEC2+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
			}
			if((idxMpT2-parent2.size())<(idxMpT1-parent1.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o1.addAll(o1.size()-3,parent2.subList(idxMpT2-2, parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC1+=Integer.parseInt(parent2.get(j).toString());					
				}
				else{
					o2.addAll(o2.size()-3,parent2.subList(idxMpT2-2,parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC2+=Integer.parseInt(parent2.get(j).toString());		
				}
			}
			if((idxMpT1-parent1.size())<(idxMpT2-parent2.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o2.addAll(o2.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC2+=Integer.parseInt(parent1.get(j).toString());	
				}
				else{
					o1.addAll(o1.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC1+=Integer.parseInt(parent1.get(j).toString());	
				}
			}
			for(int j=idxMp2-5; j>=3; j-=3){
				o2.addAll(0, parent2.subList(j,j+3));
				idxEC2+=Integer.parseInt(parent2.get(j+1).toString());
				aC2+=Integer.parseInt(parent2.get(j+1).toString());
			}
			idxEC1+=Integer.parseInt(o1.get(1).toString()) - 1;
			idxEC2+=Integer.parseInt(parent2.get(idxMp2-1).toString()) - 1;
			//Add indexes
			o1.add(0,0);
			o1.add(0,idxEC1+this.idxMS);
			o1.add(0,this.idxMS);
			o2.add(0,Integer.parseInt(parent2.get(2).toString()));
			o2.add(0,(idxEC2+Integer.parseInt(parent2.get(0).toString())));
			o2.add(0,Integer.parseInt(parent2.get(0).toString()));
			//Repair the genotypes
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP5 DONE!
		}
		else if(parent1.get(idxMp1).equals(")") && parent2.get(idxMp2).equals(":")){
			//System.out.println("OP6");
			o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
			o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
			//Crossover on the right side
			while(idxMpT1 <= parent1.size()-4 && idxMpT2 <= parent2.size()-4){
				if(random.nextBoolean()){
					o1.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					idxEC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					idxEC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					idxEC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					idxEC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1+=3;
				idxMpT2+=3;
			}
			//CrossOver the rest of configurations
			
			//Align the terminal, then crossover
			if(random.nextBoolean()){
				o1.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				o2.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				idxEC1+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
				idxEC2+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
			}
			else{
				o1.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				o2.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				idxEC1+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
				idxEC2+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
			}
			if((idxMpT2-parent2.size())<(idxMpT1-parent1.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o1.addAll(o1.size()-3,parent2.subList(idxMpT2-2, parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC1+=Integer.parseInt(parent2.get(j).toString());					
				}
				else{
					o2.addAll(o2.size()-3,parent2.subList(idxMpT2-2,parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC2+=Integer.parseInt(parent2.get(j).toString());		
				}
			}
			if((idxMpT1-parent1.size())<(idxMpT2-parent2.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o2.addAll(o2.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC2+=Integer.parseInt(parent1.get(j).toString());	
				}
				else{
					o1.addAll(o1.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC1+=Integer.parseInt(parent1.get(j).toString());	
				}
			}
			//TODO: ADD the rest of configuration to child 1
			for(int j=idxMp1-5; j>=3; j-=3){
				o1.addAll(0, parent1.subList(j,j+3));
				idxEC1+=Integer.parseInt(parent1.get(j+1).toString());
				aC1+=Integer.parseInt(parent1.get(j+1).toString());
			}
			idxEC2+=Integer.parseInt(o2.get(1).toString()) - 1;
			idxEC1+=Integer.parseInt(parent1.get(idxMp1-1).toString()) - 1;
			//Add indexes
			o1.add(0,Integer.parseInt(parent1.get(2).toString()));
			o1.add(0,(idxEC1+Integer.parseInt(parent1.get(0).toString())));
			o1.add(0,Integer.parseInt(parent1.get(0).toString()));
			o2.add(0,0);
			o2.add(0,idxEC2+this.idxMS);
			o2.add(0,this.idxMS);
			//Repair the genotypes
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP6 OK!
		}
		else if(parent1.get(idxMp1).equals(";") && parent2.get(idxMp2).equals(")")){
			//System.out.println("OP7");
			o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
			o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
			idxMpT1-=6;
			idxMpT2-=6;
			//Next crossover on the left side of complexes
			while(idxMpT1>5 && idxMpT2>5){
				if(random.nextBoolean()){
					o1.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					aC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					aC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					aC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					aC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1-=3;
				idxMpT2-=3;
			}
			//CrossOver the rest left configurations
			if(random.nextBoolean()){
				//CrossOver terminals Blocks
				for(int j=idxMpT1-2; j>=3; j-=3){
					o2.addAll(0,parent1.subList(j, j+3));
					aC2+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j>=3; j-=3){
					o1.addAll(0,parent2.subList(j,j+3));
					aC1+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			else{
				for(int j=idxMpT1-2; j>=3; j-=3){
					o1.addAll(0,parent1.subList(j, j+3));
					aC1+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j>=3; j-=3){
					o2.addAll(0,parent2.subList(j, j+3));
					aC2+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			//Add to child 2 the rest rigth configuration of parent 2
			for(int j=idxMp2+1; j<parent2.size(); j+=3){
				o2.addAll(parent2.subList(j, j+3));
				idxEC2+=Integer.parseInt(parent2.get(j+1).toString());
			}
			idxEC1+=aC1+Integer.parseInt(o1.get(o1.size()-2).toString());
			idxEC2+=aC2+Integer.parseInt(parent2.get(idxMp2-1).toString()) - 1;
			//Check if the genotype of child 1 is correct
			idxSC1=this.idxME-aC1-Integer.parseInt(o1.get(o1.size()-2).toString()) + 1;
			if(idxSC1<0){
				int restToDelete=Math.abs(idxSC1);
				int newNucleo=0;
				int j=1;
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
			//Check if the genotype of child 2 is correct
			idxSC2=this.idxMS-aC2;
			if(idxSC2<0){
				int restToDelete=Math.abs(idxSC2);
				int newNucleo=0;
				int j=1;
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
			//aC1=this.idxME-Integer.parseInt(o1.get(o1.size()-2).toString());
			//Add indexes
			o1.add(0,aC1);
			o1.add(0,this.idxME);
			o1.add(0,idxSC1);
			o2.add(0,aC2);
			o2.add(0,idxEC2+idxSC2);
			o2.add(0,idxSC2);
			//Repair the genotypes
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP7 DONE!
		}
		else if(parent1.get(idxMp1).equals(")") && parent2.get(idxMp2).equals(";")){
			//System.out.println("OP8");
			o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
			o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
			idxMpT1-=6;
			idxMpT2-=6;
			//Next crossover on the left side of complexes
			while(idxMpT1>5 && idxMpT2>5){
				if(random.nextBoolean()){
					o1.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					aC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					aC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(0,parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(0,parent2.subList(idxMpT2-2, idxMpT2+1));
					aC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					aC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1-=3;
				idxMpT2-=3;
			}
			//CrossOver the rest left configurations
			if(random.nextBoolean()){
				//Cross terminals blocks
				for(int j=idxMpT1-2; j>=3; j-=3){
					o2.addAll(0,parent1.subList(j, j+3));
					aC2+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j>=3; j-=3){
					o1.addAll(0,parent2.subList(j,j+3));
					aC1+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			else{
				for(int j=idxMpT1-2; j>=3; j-=3){
					o1.addAll(0,parent1.subList(j, j+3));
					aC1+=Integer.parseInt(parent1.get(j+1).toString());
				}
				for(int j=idxMpT2-2; j>=3; j-=3){
					o2.addAll(0,parent2.subList(j, j+3));
					aC2+=Integer.parseInt(parent2.get(j+1).toString());
				}
			}
			//Add to child 1 the rest right configuration of parent 2
			for(int j=idxMp1+1; j<parent1.size(); j+=3){
				o1.addAll(parent1.subList(j, j+3));
				idxEC1+=Integer.parseInt(parent1.get(j+1).toString());
			}
			idxEC1+=aC1+Integer.parseInt(parent1.get(idxMp1-1).toString()) - 1;
			//idxEC2+=aC2+Integer.parseInt(o2.get(o2.size()-2).toString());
			//Check if the genotype of child 1 is correct
			idxSC1=this.idxMS-aC1;
			if(idxSC1<0){
				int restToDelete=Math.abs(idxSC1);
				int newNucleo=0;
				int j=1;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o1.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						aC1-=Math.abs(lastNucleo);
						idxSC1+=Math.abs(lastNucleo);
						o1.subList(j-1,j+2).clear();
					}
					else{
						o1.set(j, newNucleo);
						aC1-=Math.abs(restToDelete);
						idxSC1+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			//Check if the genotype of child 2 is correct
			idxSC2=this.idxME-aC2-Integer.parseInt(o2.get(o2.size()-2).toString()) + 1;
			if(idxSC2<0){
				int restToDelete=Math.abs(idxSC2);
				int newNucleo=0;
				int j=1;
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
			//aC2=this.idxME-Integer.parseInt(o2.get(o2.size()-2).toString());
			//Add indexes
			o1.add(0,aC1);
			o1.add(0,idxEC1+idxSC1);
			o1.add(0,idxSC1);
			o2.add(0,aC2);
			o2.add(0,this.idxME);
			o2.add(0,idxSC2);
			//Repair genotypes
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//TODO: OP8 DONE!
		}
		else if(parent1.get(idxMp1).equals(")") && parent2.get(idxMp2).equals(")")){
			//Begin the uniform crossover on mutation origin
			if(random.nextBoolean()){
				o1.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
				o2.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
				idxMp1-=3;
				idxMp2-=3;
			}
			else{
				o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
				o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
				idxMp1-=3;
				idxMp2-=3;
			}
			idxEC1+=Integer.parseInt(o1.get(1).toString()) - 1;
			idxEC2+=Integer.parseInt(o2.get(1).toString()) - 1;
			idxMc1+=1;
			idxMc2+=1;
			//Next crossover on the left side of complexes
			while(idxMp1>5 && idxMp2>5){
				if(random.nextBoolean()){
					o1.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
					o2.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
					aC1+=Integer.parseInt(parent2.get(idxMp2-1).toString());
					aC2+=Integer.parseInt(parent1.get(idxMp1-1).toString());
					
				}
				else{
					o1.addAll(0,parent1.subList(idxMp1-2, idxMp1+1));
					o2.addAll(0,parent2.subList(idxMp2-2, idxMp2+1));
					aC1+=Integer.parseInt(parent1.get(idxMp1-1).toString());
					aC2+=Integer.parseInt(parent2.get(idxMp2-1).toString());
				}
				idxMp1-=3;
				idxMp2-=3;
				idxMc1+=3;
				idxMc2+=3;
			}
			//Align the initiators, then crossover
			if(random.nextBoolean()){
				o1.addAll(0,parent2.subList(3, 6));
				o2.addAll(0,parent1.subList(3, 6));
				aC1+=Integer.parseInt(parent2.get(4).toString());
				aC2+=Integer.parseInt(parent1.get(4).toString());
			}
			else{
				o1.addAll(0,parent1.subList(3, 6));
				o2.addAll(0,parent2.subList(3, 6));
				aC1+=Integer.parseInt(parent1.get(4).toString());
				aC2+=Integer.parseInt(parent2.get(4).toString());
			}
			idxMc1+=3;
			idxMc2+=3;
			if(idxMp2>idxMp1){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o1.addAll(3,parent2.subList(6, idxMp2+1));
					for(int j=7;j<idxMp2;j+=3){
						aC1+=Integer.parseInt(parent2.get(j).toString());
						idxMc1+=3;
					}
				}
				else{
					o2.addAll(3,parent2.subList(6, idxMp2+1));
					for(int j=7;j<idxMp2;j+=3){
						aC2+=Integer.parseInt(parent2.get(j).toString());
						idxMc2+=3;
					}
				}
			}
			if(idxMp1>idxMp2){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o2.addAll(3,parent1.subList(6, idxMp1+1));
					for(int j=7;j<idxMp1;j+=3){
						aC2+=Integer.parseInt(parent1.get(j).toString());
						idxMc2+=3;
					}
				}
				else{
					o1.addAll(3,parent1.subList(6, idxMp1+1));
					for(int j=7;j<idxMp1;j+=3){
						aC1+=Integer.parseInt(parent1.get(j).toString());
						idxMc1+=3;
					}
				}
			}
			idxEC1+=aC1;
			idxEC2+=aC2;
			//Crossover on the right side
			while(idxMpT1 <= parent1.size()-4 && idxMpT2 <= parent2.size()-4){
				if(random.nextBoolean()){
					o1.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					o2.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					idxEC1+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
					idxEC2+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					
				}
				else{
					o1.addAll(parent1.subList(idxMpT1-2, idxMpT1+1));
					o2.addAll(parent2.subList(idxMpT2-2, idxMpT2+1));
					idxEC1+=Integer.parseInt(parent1.get(idxMpT1-1).toString());
					idxEC2+=Integer.parseInt(parent2.get(idxMpT2-1).toString());
				}
				idxMpT1+=3;
				idxMpT2+=3;
			}
			//Align the terminal, then crossover
			if(random.nextBoolean()){
				o1.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				o2.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				idxEC1+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
				idxEC2+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
			}
			else{
				o1.addAll(parent1.subList(parent1.size()-3, parent1.size()));
				o2.addAll(parent2.subList(parent2.size()-3, parent2.size()));
				idxEC1+=Integer.parseInt(parent1.get(parent1.size()-2).toString());
				idxEC2+=Integer.parseInt(parent2.get(parent2.size()-2).toString());
			}
			if((idxMpT2-parent2.size())<(idxMpT1-parent1.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o1.addAll(o1.size()-3,parent2.subList(idxMpT2-2, parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC1+=Integer.parseInt(parent2.get(j).toString());					
				}
				else{
					o2.addAll(o2.size()-3,parent2.subList(idxMpT2-2,parent2.size()-3));
					for(int j=parent2.size()-5; j>idxMpT2-3; j-=3)
						idxEC2+=Integer.parseInt(parent2.get(j).toString());		
				}
			}
			if((idxMpT1-parent1.size())<(idxMpT2-parent2.size())){
				//Randomly choose to cross the remaining configuration
				if(random.nextBoolean()){
					o2.addAll(o2.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC2+=Integer.parseInt(parent1.get(j).toString());	
				}
				else{
					o1.addAll(o1.size()-3,parent1.subList(idxMpT1-2, parent1.size()-3));
					for(int j=parent1.size()-5; j>idxMpT1-3; j-=3)
						idxEC1+=Integer.parseInt(parent1.get(j).toString());	
				}
			}
			//Recalculate Initial and End Indexes
			idxSC1=this.idxMS-aC1;
			if(idxSC1<0){
				//int addNuc=random.nextInt(Math.round((Integer.parseInt(o1.get(idxMc1).toString()))));
				//idxSC1+=addNuc;
				int restToDelete=Math.abs(idxSC1);
				int newNucleo=0;
				int j=1;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o1.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					if(newNucleo<=0){
						aC1-=Math.abs(lastNucleo);
						idxSC1+=Math.abs(lastNucleo);
						o1.subList(j-1,j+2).clear();
					}
					else{
						o1.set(j, newNucleo);
						aC1-=Math.abs(restToDelete);
						idxSC1+=Math.abs(restToDelete);
					}
					restToDelete-=lastNucleo;
				}
			}
			idxSC2=this.idxMS-aC2;
			if(idxSC2<0){
				//int addNuc=random.nextInt(Math.round((Integer.parseInt(o2.get(idxMc2).toString()))));
				//idxSC2+=addNuc;
				int restToDelete=Math.abs(idxSC2);
				int newNucleo=0;
				int j=1;
				while(restToDelete>0){
					int lastNucleo=Integer.parseInt(o2.get(j).toString());
					newNucleo=lastNucleo-restToDelete;
					idxSC2+=Math.abs(newNucleo);
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
			idxEC1+=idxSC1;
			idxEC2+=idxSC2;
			o1.add(0,aC1);
			o1.add(0,idxEC1);
			o1.add(0,idxSC1);
			o2.add(0,aC2);
			o2.add(0,idxEC2);
			o2.add(0,idxSC2);
			repairGenotypeWithIndexes(o1);
			repairGenotypeWithIndexes(o2);
			//System.out.println("OP9");
			//TODO: OP9 DONE!
		}
		
		Pair<G> offspring = new Pair<G>((G) o1, (G) o2);
		return offspring;
	}
	
	private void repairGenotype(ArrayList<Object> c1, ArrayList<Object> c2){
		//Repair the child's genotype
		if(c1.get(0).equals(";")){
			c1.set(0,":");
			c1.set(2, ":");
		}
		if(c2.get(0).equals(";")){
			c2.set(0,":");
			c2.set(2, ":");
		}
		if(!c1.get(3).equals("[") && c1.get(0).equals(":")){
			c1.set(0,"[");
			c1.set(2, "]");
		}
		if(!c2.get(3).equals("[") && c2.get(0).equals(":")){
			c2.set(0,"[");
			c2.set(2, "]");
		}
		if(c1.get(3).equals("[") && c1.get(0).equals("[")){
			int tempNuc=Integer.parseInt(c1.get(1).toString());
			c1.set(4, (Integer.parseInt(c1.get(4).toString())+tempNuc));
			c1.subList(0,3).clear();
		}
		if(c2.get(3).equals("[") && c2.get(0).equals("[")){
			int tempNuc=Integer.parseInt(c2.get(1).toString());
			c2.set(4, (Integer.parseInt(c2.get(4).toString())+tempNuc));
			c2.subList(0,3).clear();
		}
		if(!c1.get(0).equals("[") && c1.get(3).equals("[")){
			c1.set(0,":");
			c1.set(2, ":");
		}
		if(!c2.get(0).equals("[") && c2.get(3).equals("[")){
			c2.set(0,":");
			c2.set(2, ":");
		}
		/*
		 * [67, 317, 17, :, 4, :, [, 13, ], (, 83, ), [, 106, ], (, 6, ), [, 13, ], (, 7, ), [, 5, ], {, 8, }, [, 2, ], <, 2, >, [, 1, ]]
		 * [62, 422, 29, [, 5, ], {, 18, }, [, 6, ], (, 144, ), [, 134, ], ;, 53, ;]
		 * Child's:
		 * [:, 4, :, [, 13, ], (, 83, ), [, 106, ], [, 1, ]]
		 * [[, 5, ], {, 18, }, [, 6, ], (, 144, ), [, 134, ], (, 6, ), [, 13, ], (, 7, ), [, 5, ], {, 8, }, [, 2, ], <, 2, >, ;, 53, ;]
		 */
		if(c1.get(c1.size()-1).equals(":") ){
			c1.set(c1.size()-3,";");
			c1.set(c1.size()-1, ";");
		}
		if(c2.get(c2.size()-1).equals(":") ){
			c2.set(c2.size()-3,";");
			c2.set(c2.size()-1, ";");
		}
		if(!c1.get(c1.size()-4).equals("]") && c1.get(c1.size()-1).equals(";")){
			c1.set(c1.size()-3,"[");
			c1.set(c1.size()-1, "]");
		}
		if(!c2.get(c2.size()-4).equals("]") && c2.get(c2.size()-1).equals(";")){
			c2.set(c2.size()-3,"[");
			c2.set(c2.size()-1, "]");
		}
		if(c1.get(c1.size()-4).equals("]") && c1.get(c1.size()-1).equals("]")){
			int tempNuc=Integer.parseInt(c1.get(c1.size()-2).toString());
			c1.set(c1.size()-5, (Integer.parseInt(c1.get(c1.size()-5).toString())+tempNuc));
			c1.subList(c1.size()-3, c1.size()).clear();
			
		}
		if(c2.get(c2.size()-4).equals("]") && c2.get(c2.size()-1).equals("]")){
			int tempNuc=Integer.parseInt(c2.get(c2.size()-2).toString());
			c2.set(c2.size()-5, (Integer.parseInt(c2.get(c2.size()-5).toString())+tempNuc));
			c2.subList(c2.size()-3, c2.size()).clear();
		}
	}
	
	private void repairGenotypeWithIndexes(ArrayList<Object> c){
		if(c.get(3).equals(";")){
			c.set(3,":");
			c.set(5, ":");
		}
		if(!c.get(3).equals("[") && !c.get(3).equals("[") && 
				Integer.parseInt(c.get(2).toString()) != 0){
			c.set(3,"[");
			c.set(5, "]");
		}
		if(!c.get(6).equals("[") && c.get(3).equals(":") && 
				Integer.parseInt(c.get(2).toString()) == 0){
			c.set(6,"[");
			c.set(8, "]");
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
		if(!c.get(c.size()-4).equals("]") && !c.get(c.size()-1).equals("]") &&
				Integer.parseInt(c.get(1).toString()) != this.idxME){
			c.set(c.size()-3,"[");
			c.set(c.size()-1, "]");
		}
		if(!c.get(c.size()-4).equals("]") && c.get(c.size()-1).equals(";") &&
				Integer.parseInt(c.get(1).toString()) == this.idxME){
			c.set(c.size()-6,"[");
			c.set(c.size()-4, "]");
		}
		if(c.get(c.size()-4).equals("]") && c.get(c.size()-1).equals("]")){
			int tempNuc=Integer.parseInt(c.get(c.size()-2).toString());
			c.set(c.size()-5, (Integer.parseInt(c.get(c.size()-5).toString())+tempNuc));
			c.subList(c.size()-3, c.size()).clear();
			
		}
	}
}