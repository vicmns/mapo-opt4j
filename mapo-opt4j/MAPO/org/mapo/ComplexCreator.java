package org.mapo;
import java.util.ArrayList;
import java.util.Random;

import org.opt4j.core.problem.Creator;
import org.opt4j.genotype.DynamicListGenotype;

import com.google.inject.Inject;


public class ComplexCreator implements Creator<DynamicListGenotype<Object>>{
	
	protected final ComplexProblem problem;
	protected final Random rand;
	
	@Inject
	public ComplexCreator(ComplexProblem problem, Random rand){
		this.problem=problem;
		this.rand=rand;
	}
		
	private String getCloseChar(char openChar){
		if(openChar=='[')
			return "]";
	 	if(openChar=='(')
	 		return ")";
	 	if(openChar=='<')
	 		return ">";
	 	if(openChar=='{')
	 		return "}";
	 	if(openChar==':')
	 		return ":";
	 	if(openChar==';')
	 		return ";";
	 	else
	 		return "";
	}
	
	public ArrayList<Object> complexGenerator(int idxMS, int idxME){
		ArrayList<Object> complex = new ArrayList<Object>();
		ArrayList<Object> tempList = new ArrayList<Object>();
		int centerNuc=0;
		int nNucleo=0;
		int len=0;
		int idxS=0; //Start nucleotide in DNA/cDNA
		int idxE=0; //End nucleotide in DNA/cDNA
		int nNucAvLeft=0;
		int nNucAvRigth=0;
		int alpha=0;
		int length=problem.getABLength();
		//First tuple {i,j,alpha} of individual representing the substring
		do{
			idxS = rand.nextInt(idxMS);
			idxE = rand.nextInt(length);
			len=idxE-idxS;
		}while(idxS>=idxMS  || idxME>=idxE);
		complex.add(idxS);
		complex.add(idxE);
		/*Next tuple representing the configuration
		  Initial configuration */
		if(rand.nextBoolean()){
			nNucleo=0;
			do{
				nNucleo=rand.nextInt((len/2)+2);
			}while(nNucleo < 3);
			nNucAvLeft=Math.abs(idxMS-idxS-rand.nextInt(nNucleo));
			nNucAvRigth=len-(nNucAvLeft+nNucleo);
			centerNuc=nNucleo;
		}
		else{
			nNucleo = idxME - idxMS + 1;
			nNucAvLeft=idxMS-idxS;
			nNucAvRigth=idxE-idxME-1;
			centerNuc=nNucleo;
		}
		int lastNNucAvL=nNucAvLeft;
		/*Interior representation filled
		 	from center to left*/
		while(nNucAvLeft>0){
			//Add complementary section after a random non complementary configuration
			do{
				nNucleo=rand.nextInt(nNucAvLeft+1);
			}while(nNucleo<1);
			tempList.add(0,"]");
			tempList.add(0,nNucleo);
			tempList.add(0,"[");
			nNucAvLeft-=nNucleo;
			//Random non-complementary configuration selection
			if(nNucAvLeft>0){
				char[] opLenguage = {':','(','<','{'};
				char op;
				do{
					op = opLenguage[rand.nextInt(4)];
				}while(op=='(' && nNucAvLeft < 3);
				switch(op){
				case ':':
					do{
						nNucleo = rand.nextInt(nNucAvLeft+1);
					}while(nNucleo<1);
					/*
					 * Check if there are any nucelotides aviable before adding the overhang
					 * if so, add the rest nucloe avaible to the last configuration 
					 */
					if((nNucAvLeft-nNucleo) > 0){
						int lastNuc=Integer.parseInt(tempList.get(1).toString());
						tempList.set(1, lastNuc+(nNucAvLeft-nNucleo));
						nNucAvLeft-=(nNucAvLeft-nNucleo);
					}
					break;
				case '(':
					if(nNucleo>3)
						do{
							nNucleo = rand.nextInt(nNucAvLeft+1);
						}while(nNucleo<3);
					break;
				case '<':
					do{
						nNucleo = rand.nextInt(nNucAvLeft+1);
					}while(nNucleo<1);
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvLeft+1);
					}while(nNucleo<1);
					break;
				}
				tempList.add(0,getCloseChar(op));
				tempList.add(0,nNucleo);
				tempList.add(0,""+op);
				nNucAvLeft-=nNucleo;
			}
		}
		/*
		 * Check if last configuration is valid (":#:" or "[#]")
		 */
		if(!tempList.isEmpty()){
			int lastElement=0;
			if(!tempList.get(lastElement).equals(":") && !tempList.get(lastElement).equals("[")){
				int lastNuc=Integer.parseInt(tempList.get(lastElement+1).toString());
				tempList.remove(lastElement);
				tempList.remove(lastElement);
				tempList.remove(lastElement);
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					tempList.add(0,":");
					tempList.add(0,lastNuc);
					tempList.add(0,":");
				}
				else{
					//Complementary Section
					if(tempList.size()>1){
						int newNuc=Integer.parseInt(tempList.get(1).toString());
						newNuc+=lastNuc;
						tempList.set(1, newNuc);
					}
					else{
						tempList.add(0,"]");
						tempList.add(0,lastNuc);
						tempList.add(0,"[");
					}
				}
			}
		}
		alpha=lastNNucAvL;
		complex.add(alpha);
		complex.addAll(tempList);
		tempList=new ArrayList<Object>();
		/*Interior representation filled
	 		from center to rigth*/
		while(nNucAvRigth>0){
			//Add complementary section after a random non complementary configuration
			do{
				nNucleo=rand.nextInt(nNucAvRigth+1);
			}while(nNucleo<1);
			tempList.add("[");
			tempList.add(nNucleo);
			tempList.add("]");
			nNucAvRigth-=nNucleo;
			//Random non-complementary configuration selection
			if(nNucAvRigth>0){
				char[] opLenguage = {';','(','<','{'};
				char op;
				do{
					op = opLenguage[rand.nextInt(4)];
				}while(op=='(' && nNucAvRigth < 3);
				switch(op){
				case ';':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					/*
					 * Check if there are any nucelotides aviable before adding the overhang
					 * if so, add the rest nucleo avaible to the last configuration 
					 */
					int diffNuc=nNucAvRigth-nNucleo;
					if(diffNuc > 0){
						int lastNuc=Integer.parseInt(tempList.get(tempList.size()-2).toString());
						int sumNuc= lastNuc+diffNuc;
						tempList.set(tempList.size()-2, sumNuc);
						nNucAvRigth-=diffNuc;
					}
					break;
				case '(':
					if(nNucleo>3)
						do{
							nNucleo = rand.nextInt(nNucAvRigth+1);
						}while(nNucleo<3);
					break;
				case '<':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					break;
				}
				tempList.add(""+op);
				tempList.add(nNucleo);
				tempList.add(getCloseChar(op));
				nNucAvRigth-=nNucleo;
			}
		}
		/*
		 * Check if last configuration is valid (";#;" or "[#]")
		 */
		if(!tempList.isEmpty()){
			int lastElement=tempList.size()-1;
			if(!tempList.get(lastElement).equals(";") || !tempList.get(lastElement).equals("]")){
				int lastNuc=Integer.parseInt(tempList.get(lastElement-1).toString());
				tempList.remove(lastElement);
				tempList.remove(lastElement-1);
				tempList.remove(lastElement-2);
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					tempList.add(";");
					tempList.add(lastNuc);
					tempList.add(";");
				}
				else{
					//Complementary Section
					if(tempList.size()>1){
						int newNuc=Integer.parseInt(tempList.get(tempList.size()-2).toString());
						newNuc+=lastNuc;
						tempList.set(tempList.size()-2, newNuc);
					}
					else{
						tempList.add("[");
						tempList.add(lastNuc);
						tempList.add("]");
					}
				}
			}
		}
		complex.add("(");
		complex.add(centerNuc);
		complex.add(")");
		complex.addAll(tempList);
		return complex;
	}
	
	public DynamicListGenotype<Object> create(){
		DynamicListGenotype<Object> genotype = new DynamicListGenotype<Object>();
		genotype.addAll(complexGenerator(problem.getIdxS(),problem.getIdxE()));
		return genotype;
	}
}
