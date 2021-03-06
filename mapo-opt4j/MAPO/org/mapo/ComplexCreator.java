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
		switch(rand.nextInt(3)){
		case 0: //Interior loop representation
			//System.out.println("Generating individual at Case 0");
			return interiorLoopMutation(idxMS, idxME);
		case 1:
			//System.out.println("Generating individual at Case 1");
			return leftOverhangMutation(idxMS, idxME);
		case 2:
			//System.out.println("Generating individual at Case 2");
			return rightOverhangMutation(idxMS, idxME);
		default:
			//System.out.println("Generating individual at Case 0");
			return interiorLoopMutation(idxMS, idxME);
		}
	}
	
	public DynamicListGenotype<Object> create(){
		DynamicListGenotype<Object> genotype = new DynamicListGenotype<Object>();
		genotype.addAll(complexGenerator(problem.getIdxS(),problem.getIdxE()));
		return genotype;
	}
	
	private ArrayList<Object> interiorLoopMutation(int idxMS, int idxME){
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
		int numOverhangs = 0;
		//int numInfBulgue = 0;
		//int numSupBulgue = 0;
		boolean hasInfBulgue = false;
		//First tuple {i,j,alpha} of individual representing the substring
		do{
			idxS = rand.nextInt(idxMS);
			idxE = rand.nextInt(length);
			len = idxE - idxS + 1;
		}while(idxS >= idxMS + 3  || idxME + 3 >= idxE);
		complex.add(idxS);
		complex.add(idxE);
		/*Next tuple representing the configuration
		  Initial configuration */
		if(rand.nextBoolean()){
			nNucleo=0;
			do{
				nNucleo=rand.nextInt((len))-(idxME-idxMS + 1);
			}while(nNucleo < (idxME-idxMS+1) || nNucleo > (len-(idxME-idxMS+1)));
			do{
				nNucAvLeft=Math.abs(idxMS-idxS-rand.nextInt(nNucleo-1)) + 1;
				nNucAvRigth=len-(nNucAvLeft+nNucleo);
			}while(nNucAvLeft<1 || nNucAvRigth<1);
			centerNuc=nNucleo;
		}
		else{
			nNucleo = idxME - idxMS + 1;
			nNucAvLeft = idxMS - idxS;
			nNucAvRigth = idxE - idxME;
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
					 * if so, add the rest nucleo aviable to the last configuration 
					 */
					if((nNucAvLeft-nNucleo) > 0){
						int lastNuc=Integer.parseInt(tempList.get(1).toString());
						tempList.set(1, lastNuc+(nNucAvLeft-nNucleo));
						nNucAvLeft-=(nNucAvLeft-nNucleo);
					}
					numOverhangs+=nNucleo;
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
					//numSupBulgue += nNucleo;
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvLeft+1);
					}while(nNucleo<1);
					//numInfBulgue+=nNucleo;
					hasInfBulgue = true;
					break;
				}
				tempList.add(0,getCloseChar(op));
				tempList.add(0,nNucleo);
				tempList.add(0,""+op);
				if(!hasInfBulgue)
					nNucAvLeft-=nNucleo;
				hasInfBulgue = false;
			}
		}
		/*
		 * Check if last configuration is valid (":#:" or "[#]")
		 */
		if(!tempList.isEmpty()){
			int lastElement=0;
			if(!tempList.get(lastElement).equals(":") && !tempList.get(lastElement).equals("[")){
				int lastNuc=Integer.parseInt(tempList.get(lastElement+1).toString());
				if(tempList.get(lastElement).equals("{")){
					//numInfBulgue -= lastNuc;
				}
				if(tempList.get(lastElement).equals("<")){
					//numSupBulgue -= lastNuc;
				}
				tempList.remove(lastElement);
				tempList.remove(lastElement);
				tempList.remove(lastElement);
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					tempList.add(0,":");
					tempList.add(0,lastNuc);
					tempList.add(0,":");
					numOverhangs += lastNuc;
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
		//nNucAvLeft = 0;
		/*if(numInfBulgue > 0){
			nNucAvLeft = numInfBulgue;
		}
		else
			nNucAvLeft = 0*/
		alpha=lastNNucAvL - nNucAvLeft;
		complex.add(alpha);
		complex.set(0, idxS + nNucAvLeft );
		complex.addAll(tempList);
		tempList=new ArrayList<Object>();
		/*Interior representation filled
	 		from center to rigth*/
		//nNucAvRigth += (nNucAvLeft - numSupBulgue);
		//numInfBulgue = 0;
		//numSupBulgue = 0;
		//hasInfBulgue = false;
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
					numOverhangs+=nNucleo;
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
					//numSupBulgue += nNucleo;
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					//numInfBulgue+=nNucleo;
					hasInfBulgue = true;
					break;
				}
				tempList.add(""+op);
				tempList.add(nNucleo);
				tempList.add(getCloseChar(op));
				if(!hasInfBulgue)
					nNucAvRigth-=nNucleo;
				hasInfBulgue = false;
			}
		}
		/*
		 * Check if last configuration is valid (";#;" or "[#]")
		 */
		if(!tempList.isEmpty()){
			int lastElement=tempList.size()-1;
			if(!tempList.get(lastElement).equals(";") && !tempList.get(lastElement).equals("]")){
				int lastNuc=Integer.parseInt(tempList.get(lastElement-1).toString());
				if(tempList.get(lastElement).equals("}")){
					//numInfBulgue -= lastNuc;
				}
				if(tempList.get(lastElement).equals(">")){
					//numSupBulgue -= lastNuc;
				}
				tempList.subList(lastElement-2,lastElement+1).clear();
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					tempList.add(";");
					tempList.add(lastNuc);
					tempList.add(";");
					numOverhangs += lastNuc;
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
		nNucAvRigth=0;
		/*if(numInfBulgue > 0 && ( (numOverhangs + numSupBulgue) > numInfBulgue)){
			//nNucAvRigth +=  numInfBulgue;
			complex.set(1, Integer.parseInt(complex.get(1).toString()) - numInfBulgue);
		}
		else if(numInfBulgue > 0 && ( (numOverhangs + numSupBulgue) < numInfBulgue)){
			nNucAvRigth = numInfBulgue - (numOverhangs + numSupBulgue);
			//complex.set(1, Integer.parseInt(complex.get(1).toString()) - numInfBulgue - numSupBulgue + nNucAvRigth);
			complex.set(1, Integer.parseInt(complex.get(1).toString()) - numInfBulgue);
		}*/		
		complex.add("(");
		complex.add(centerNuc);
		complex.add(")");
		complex.addAll(tempList);
		tempList.clear();
		return complex;
	}
	
	private ArrayList<Object> leftOverhangMutation(int idxMS, int idxME){
 		ArrayList<Object> complex = new ArrayList<Object>();
		ArrayList<Object> tempList = new ArrayList<Object>();
		int nNucleo=0;
		int len=0;
		int idxS = idxMS; //Start nucleotide in DNA/cDNA
		int idxE=0; //End nucleotide in DNA/cDNA
		int nNucAvRigth=0;
		int length=problem.getABLength();
		int numOverhangs = 0;
		int numInfBulgue = 0;
		int numSupBulgue = 0;
		boolean hasInfBulgue = false;
		//First tuple {i,j,alpha} of individual representing the substring
		do{
			idxE = rand.nextInt(length);
			len = idxE - idxS + 1;
		}while(idxME + 3 >= idxE);
		complex.add(idxS);
		complex.add(idxE);
		complex.add(0); //Alpha position
		int leftOverhang=0;
		do{
			leftOverhang= rand.nextInt(8) + 1;
		}while(leftOverhang < idxME-idxMS + 1 || leftOverhang>=len);
		complex.add(":"); complex.add(leftOverhang); complex.add(":");
		numOverhangs+=leftOverhang;
		nNucAvRigth = len - leftOverhang;
		while(nNucAvRigth>0){
			//Add complementary section after a random non complementary configuration
			do{
				nNucleo=rand.nextInt(nNucAvRigth + 1);
			}while(nNucleo < 1);
			tempList.add("[");
			tempList.add(nNucleo);
			tempList.add("]");
			nNucAvRigth -= nNucleo;
			//Random non-complementary configuration selection
			if(nNucAvRigth > 0){
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
					numOverhangs+=nNucleo;
					break;
				case '(':
						do{
							nNucleo = rand.nextInt(nNucAvRigth+1);
						}while(nNucleo<3);
					break;
				case '<':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					numSupBulgue += nNucleo;
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					numInfBulgue+=nNucleo;
					hasInfBulgue = true;
					//nNucleo=0;
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
			if(!tempList.get(lastElement).equals(";") && !tempList.get(lastElement).equals("]")){
				int lastNuc=Integer.parseInt(tempList.get(lastElement-1).toString());
				if(tempList.get(lastElement).equals("}"))
					hasInfBulgue=false;
				if(tempList.get(lastElement).equals(">"))
					numSupBulgue -= lastNuc;
				tempList.subList(lastElement-2, lastElement+1).clear();
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					tempList.add(";");
					tempList.add(lastNuc);
					tempList.add(";");
					numOverhangs += lastNuc;
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
		nNucAvRigth=0;
		if(hasInfBulgue && ( (numOverhangs + numSupBulgue) > numInfBulgue)){
			//nNucAvRigth +=  numInfBulgue;
			complex.set(1, Integer.parseInt(complex.get(1).toString()) - numInfBulgue);
		}
		else if(hasInfBulgue && ( (numOverhangs + numSupBulgue) < numInfBulgue)){
			nNucAvRigth = numInfBulgue - (numOverhangs + numSupBulgue);
			complex.set(1, Integer.parseInt(complex.get(1).toString())-numInfBulgue);
			//complex.set(1, Integer.parseInt(complex.get(1).toString()) - numInfBulgue  + nNucAvRigth);
																		// - numInfBulgue - numOverhangs - numSupBulgue 
		}		
		complex.addAll(tempList);
		tempList.clear();
		return complex;
	}

	private ArrayList<Object> rightOverhangMutation(int idxMS, int idxME){
		ArrayList<Object> complex = new ArrayList<Object>();
		ArrayList<Object> tempList = new ArrayList<Object>();
		int nNucleo=0;
		int len=0;
		int idxS=0; //Start nucleotide in DNA/cDNA
		int idxE=idxME; //End nucleotide in DNA/cDNA 
		int nNucAvLeft=0;
		int numOverhangs=0;
		int numInfBulgue=0;
		int numSupBulgue=0;
		boolean hasInfBulgue = false;
		//First tuple {i,j,alpha} of individual representing the substring
		do{
			idxS = rand.nextInt(idxMS);
			len=idxE-idxS+1;
		}while(idxS>=idxMS-3);
		//Generate right overhang
		int rightOverhang=0;
		do{
			rightOverhang=rand.nextInt(8)+1;
		}while(rightOverhang<(idxME-idxMS+1) || rightOverhang>=len);
		complex.add(";"); complex.add(rightOverhang); complex.add(";");
		numOverhangs+=rightOverhang;
		int alpha=nNucAvLeft=len-rightOverhang;
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
					numOverhangs+=nNucleo;
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
					numSupBulgue += nNucleo;
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvLeft+1);
					}while(nNucleo<1);
					numInfBulgue+=nNucleo;
					hasInfBulgue=true;
					//nNucleo=0;
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
				if(tempList.get(lastElement).equals("{")){
					hasInfBulgue=false;
					numInfBulgue -= lastNuc;
				}
				tempList.subList(lastElement, lastElement+3).clear();
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					tempList.add(0,":");
					tempList.add(0,lastNuc);
					tempList.add(0,":");
					numOverhangs += lastNuc;
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
		if(numInfBulgue > 0 ){
			//nNucAvLeft +=  numInfBulgue;
			complex.add(0,idxS + numInfBulgue );
			alpha -= numInfBulgue;
		}		
		else{
			complex.add(0,idxS);
			//alpha -= numInfBulgue;
		}
		complex.add(1,idxE);
		complex.add(2,alpha);
		complex.addAll(3,tempList);
		tempList.clear();
		return complex;
	}

}