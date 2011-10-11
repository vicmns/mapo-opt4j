package org.mapo;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.opt4j.core.problem.Creator;
import org.opt4j.genotype.SelectGenotype;


public class complexCreator implements Creator<SelectGenotype<String>>{
	private int idxSs;
	private int idxEe;
	private int alpha;
	private int idxMS;
	private int idxME;
	Random rand = new Random();
	private String Gene; //Correct cDNA/DNA Sequence
	private String AB; //Mutated cDNA/RNA
	
	public  int[] findMutation(String a, String b){
		int firstPos=0;
		int[] idx = new int[]{-1,-1};
		int len=0;
		boolean hasMutation=true;
		int i=0;
		if(a.length()>b.length())
			len=b.length();
		else
			len=a.length();
		while(hasMutation && i<len){
			if(a.charAt(i)!=b.charAt(i)){
				if(idx[0]<0)
					idx[0]=i;
				if(a.charAt(i+1)==b.charAt(i+1)){
					hasMutation=false;
				}
				firstPos++;
			}
			i++;
		}
		idx[1]=idx[0]+firstPos-1;
		return idx;
	}
	
	public  String sequenceAligment(String seqA, String seqB){
		String aligment="";
		int len=0;
		if(seqA.length()>=seqB.length()){
			len=seqA.length();
			for(int i=0; i<len; i++){
				if(i<seqB.length()){
					if(seqA.charAt(i)==seqB.charAt(i)){
						aligment+=seqA.charAt(i);
					}
					else
						aligment+="-";
				}
				else{
					aligment+="-";
				}
			}
		}
		else{
			len=seqB.length();
			for(int i=0; i<=len; i++){
				if(i<seqA.length()){
					if(seqA.charAt(i)==seqB.charAt(i)){
						aligment+=seqB.charAt(i);
					}
					else
						aligment+="-";
				}
				else{
					aligment+="-";
				}
			}
		}
		return aligment;
	}
	
	public SelectGenotype<String> create(){
		//StringBuffer complex = new StringBuffer();
		List<String> complex = new ArrayList<String>();
		ArrayList<Integer> lSideNuc = new ArrayList<Integer>();
		ArrayList<Integer> rSideNuc = new ArrayList<Integer>();
		ArrayList<String> lSideStruct = new ArrayList<String>();
		ArrayList<String> rSideStruct = new ArrayList<String>();
		int centerNuc=0;
		char centerStruc='(';
		int nNucleo=0;
		int len=0;
		int idxS=0; //Start nucleotide in DNA/cDNA
		int idxE=0; //End nucleotide in DNA/cDNA
		int nNucAvLeft=0;
		int nNucAvRigth=0;
		int length=AB.length();
		boolean iniFlag=true;
		String struct = "";
		//First tuple {i,j,alpha} of individual representing the substring
		do{
			idxS = rand.nextInt(idxMS);
			idxE = rand.nextInt(length);
			len=idxE-idxS;
		}while(idxS>=idxMS  || idxME>=idxE);
		idxSs=idxS;
		idxEe=idxE;
		System.out.println("Start Index: "+idxS+" End Index: "+idxE+" Length: "+len);
		complex.add(idxSs+","+idxEe);
		/*Next tuple representing the configuration
		  Initial configuration */
		if(rand.nextBoolean()){
			nNucleo=0;
			do{
				nNucleo=rand.nextInt(len/2);
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
			lSideNuc.add(nNucleo);
			lSideStruct.add("[");
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
					 * if so, add the rest nuclo avaible to the last configuration 
					 */
					if((nNucAvLeft-nNucleo) > 0){
						int lastNuc=lSideNuc.get(lSideNuc.size()-1);
						lSideNuc.set(lSideNuc.size()-1, lastNuc+(nNucAvLeft-nNucleo));
						nNucAvLeft-=(nNucAvLeft-nNucleo);
					}
					lSideNuc.add(nNucleo);
					lSideStruct.add(""+op);
					break;
				case '(':
					if(nNucleo>3)
						do{
							nNucleo = rand.nextInt(nNucAvLeft+1);
						}while(nNucleo<3);
					lSideNuc.add(nNucleo);
					lSideStruct.add(""+op);
					break;
				case '<':
					do{
						nNucleo = rand.nextInt(nNucAvLeft+1);
					}while(nNucleo<1);
					lSideNuc.add(nNucleo);
					lSideStruct.add(""+op);
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvLeft+1);
					}while(nNucleo<1);
					lSideNuc.add(nNucleo);
					lSideStruct.add(""+op);
					break;
				}
				nNucAvLeft-=nNucleo;
			}
		}
		/*
		 * Check if last configuration is valid (":#:" or "[#]")
		 */
		if(!lSideStruct.isEmpty()){
			int lastElement=lSideStruct.size()-1;
			if(!lSideStruct.get(lastElement).equals(":") && !lSideStruct.get(lastElement).equals("[")){
				int lastNuc=lSideNuc.get(lastElement);
				lSideNuc.remove(lastElement);
				lSideStruct.remove(lastElement);
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					lSideNuc.add(lastNuc);
					lSideStruct.add(":");
				}
				else{
					/*
					 * TODO repair this section, not always will be a complementary section
					 * before any other
					 */
					//Complementary Section
					if(!lSideNuc.isEmpty()){
						int newNuc=lSideNuc.get(lSideNuc.size()-1);
						newNuc+=lastNuc;
						lSideNuc.set(lSideNuc.size()-1, newNuc);
					}
					else{
						lSideNuc.add(lastNuc);
						lSideStruct.add("[");
					}
				}
			}
		}
		this.alpha=idxSs+(lastNNucAvL);
		complex.add(","+this.alpha);
		/*Interior representation filled
	 		from center to rigth*/
		while(nNucAvRigth>0){
			//Add complementary section after a random non complementary configuration
			do{
				nNucleo=rand.nextInt(nNucAvRigth+1);
			}while(nNucleo<1);
			rSideNuc.add(nNucleo);
			rSideStruct.add("[");
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
						int lastNuc=rSideNuc.get(rSideNuc.size()-1);
						int sumNuc= lastNuc+diffNuc;
						rSideNuc.set(rSideNuc.size()-1, sumNuc);
						nNucAvRigth-=diffNuc;
					}
					rSideNuc.add(nNucleo);
					rSideStruct.add(""+op);
					break;
				case '(':
					if(nNucleo>3)
						do{
							nNucleo = rand.nextInt(nNucAvRigth+1);
						}while(nNucleo<3);
					rSideNuc.add(nNucleo);
					rSideStruct.add(""+op);
					break;
				case '<':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					rSideNuc.add(nNucleo);
					rSideStruct.add(""+op);
					break;
				case '{':
					do{
						nNucleo = rand.nextInt(nNucAvRigth+1);
					}while(nNucleo<1);
					rSideNuc.add(nNucleo);
					rSideStruct.add(""+op);
					break;
				}
				nNucAvRigth-=nNucleo;
			}
		}
		/*
		 * Check if last configuration is valid (";#;" or "[#]")
		 */
		if(!rSideStruct.isEmpty()){
			int lastElement=rSideStruct.size()-1;
			if(!rSideStruct.get(lastElement).equals(";") && !rSideStruct.get(lastElement).equals("[")){
				int lastNuc=rSideNuc.get(lastElement);
				rSideNuc.remove(lastElement);
				rSideStruct.remove(lastElement);
				//Select randomly the next correct configuration
				if(rand.nextBoolean()){
					//Left Overhang
					rSideNuc.add(lastNuc);
					rSideStruct.add(";");
				}
				else{
					//Complementary Section
					if(!lSideNuc.isEmpty()){
						int newNuc=rSideNuc.get(rSideNuc.size()-1);
						newNuc+=lastNuc;
						rSideNuc.set(rSideNuc.size()-1, newNuc);
					}
					else{
						rSideNuc.add(lastNuc);
						rSideStruct.add("[");
					}
				}
			}
		}
		/*
		 * Generate the complex String
		 */
		for(int i=lSideStruct.size()-1; i>=0; i--){
			String buffer="";
			buffer+=lSideStruct.get(i)+lSideNuc.get(i)+getCloseChar(lSideStruct.get(i));
			complex.add(buffer);
		}
		complex.add(""+centerStruc+centerNuc+getCloseChar(centerStruc));
		for(int i=0; i<rSideStruct.size(); i++){
			String buffer="";
			buffer+=rSideStruct.get(i)+rSideNuc.get(i)+getCloseChar(rSideStruct.get(i));
			complex.add(buffer);
		}
		SelectGenotype<String> genotype = new SelectGenotype<String>(complex);
		return genotype;
		//return complex.toString();
	}
	
	private String getCloseChar(String openChar){
		if(openChar.equals("["))
			return "]";
	 	if(openChar.equals("("))
	 		return ")";
	 	if(openChar.equals("<"))
	 		return ">";
	 	if(openChar.equals("{"))
	 		return "}";
	 	if(openChar.equals(":"))
	 		return ":";
	 	if(openChar.equals(";"))
	 		return ";";
	 	else
	 		return "";
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
}
