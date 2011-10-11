package org.mapo;

import java.util.ArrayList;
import java.util.Random;

import org.opt4j.common.random.Rand;
import org.opt4j.genotype.ListGenotype;
import org.opt4j.genotype.PermutationGenotype;
import org.opt4j.genotype.SelectGenotype;
import org.opt4j.operator.mutate.MutationRate;

import com.google.inject.Inject;

public class MutateInversionComplex implements MutateComplex {
	
	protected final Random random;
	
	protected final MutationRate mutationRate;
	
	private Object[] pSizeParse;
	private Object[] pSignParse;
	private int idxSp;
	private int idxEp;
	private int alphaP;
	
	@Inject
	public MutateInversionComplex(final MutationRate mutationRate, Rand random) {
		this.mutationRate = mutationRate;
		this.random = random;
	}

	@Override
	@SuppressWarnings("unchecked")
	public void mutate(ListGenotype<?> genotype, double p) {
		ListGenotype<Object> mutation = (ListGenotype<Object>) genotype;
		Object[][] pParse=complexParse(toString(genotype.toArray()));
		this.pSizeParse=pParse[0];
		this.pSignParse=pParse[1];
		this.idxSp=Integer.parseInt(pParse[2][0].toString());
		this.idxEp=Integer.parseInt(pParse[2][1].toString());
		this.alphaP=Integer.parseInt(pParse[2][2].toString());
		
		ArrayList<String> c= new ArrayList<String>();
		ArrayList<String> newSizeComplex= new ArrayList<String>();
		int fAlphaP=this.idxSp;
		int idxMP=0; //Mutation index on parent
		int idxSC=this.idxSp;
		int idxEC=this.idxEp;
		int alphaC=this.alphaP;
		int i=0;
		/*
		 * Find where the mutation starts
		 */
		while(fAlphaP!=this.alphaP){
			fAlphaP+=Integer.parseInt(pSizeParse[i].toString());
			i++;
			idxMP=i;
		}
		/*
		 * Select where to do the inversion (right or left)
		 */
		if(random.nextBoolean()){			
			/*
			 * Select a random position to the left to make the inversion
			 */
			alphaC=0;
			int invIdx=random.nextInt(idxMP);
			for(int j=0; j<invIdx; j++){
				alphaC+=Integer.parseInt(pSizeParse[j].toString());
				newSizeComplex.add(pSizeParse[j].toString());
			}
			for(int j=idxMP; j>=invIdx; j--){
				//Calculate the new alpha
				alphaC+=Integer.parseInt(pSizeParse[j].toString());
				newSizeComplex.add(pSizeParse[j].toString());
			}
			alphaC-=Integer.parseInt(pSizeParse[invIdx].toString());
			//Add the rest of the complex
			for(int j=idxMP+1; j<pSizeParse.length;j++){
				newSizeComplex.add(pSizeParse[j].toString());
			}
			for(int j=0; j<pSizeParse.length; j++){
				c.add(getOpenChar(pSignParse[j].toString())+newSizeComplex.get(j)+pSignParse[j].toString());
			}
			c.add(0,idxSC+","+idxEC+","+(alphaC+idxSC));
		}
		else{
			//Add the first part of the complex
			for(int j=0; j<idxMP;j++){
				c.add(getOpenChar(pSignParse[j].toString())+pSizeParse[j].toString()+pSignParse[j].toString());
			}
			/*
			 * Select a random position to the right to make the inversion
			 */
			int invIdx=0;
			do{
				invIdx=random.nextInt(pSizeParse.length);
			}while(invIdx<idxMP || invIdx>=pSizeParse.length);
			i=idxMP;
			for(int j=invIdx; j>=idxMP; j--){
				c.add(getOpenChar(pSignParse[i].toString())+pSizeParse[j].toString()+pSignParse[i].toString());
				i++;
			}
			for(int j=invIdx+1; j<pSizeParse.length; j++){
				c.add(getOpenChar(pSignParse[j].toString())+pSizeParse[j].toString()+pSignParse[j].toString());
			}
			c.add(0,idxSC+","+idxEC+","+(alphaC));
		}
		mutation.addAll(c);
		//return toString(c.toArray());
	}
	
	public Object[][] complexParse(String complex){
		ArrayList<String> parse = new ArrayList<String>();
		ArrayList<String> sParse = new ArrayList<String>();
		ArrayList<String> idxSE=new ArrayList<String>();
		String temp="";
		String signTemp="";
		String idxTemp="";
		boolean next=false;
		boolean nextColon=false;
		boolean nextSemicolon=false;
		boolean nextIsComplex=false;
		for(int i=0; i<complex.length(); i++){
			if(!nextIsComplex){
				if(Character.isDigit(complex.charAt(i))){
					idxTemp+=complex.charAt(i);
				}
				if(complex.charAt(i)==','){
					idxSE.add(idxTemp);
					idxTemp="";
				}
				if(complex.charAt(i+1)==':' || complex.charAt(i+1)=='['){				
					idxSE.add(idxTemp);
					nextIsComplex=true;
				}
			}
			else{
				if(complex.charAt(i)==':'){
					if(nextColon){
						next=true;
						nextColon=false;
						signTemp+=complex.charAt(i);
					}
					else
						nextColon=true;
				}
				if(Character.isDigit(complex.charAt(i))){
					temp+=complex.charAt(i);
				}
				
				if(complex.charAt(i)==']' || complex.charAt(i)==')' 
					|| complex.charAt(i)=='>' || complex.charAt(i)=='}'){
					next=true;
					signTemp+=complex.charAt(i);
				}
				if(complex.charAt(i)==';'){
					if(nextSemicolon){
						next=true;
						nextSemicolon=false;
						signTemp+=complex.charAt(i);
					}
					else
						nextSemicolon=true;
				}
				if(next){
					parse.add(temp);
					sParse.add(signTemp);
					temp="";
					signTemp="";
					next=false;
				}
			}
		}
		//System.out.println("The Nucleotides ArrayList contains: "+parse);
		//System.out.println("The Sign ArrayList contains: "+sParse);
		Object[][] complexParse= new Object[3][];
		complexParse[0]=parse.toArray();
		complexParse[1]=sParse.toArray();
		complexParse[2]=idxSE.toArray();
		return complexParse;
	}
	
	private String toString(Object[] arrayObject){
		String str="";
		for(int i=0; i<arrayObject.length; i++){
			str+=arrayObject[i];
		}
		return str;
	}
	
	private String getOpenChar(String openChar){
		if(openChar.equals("]"))
			return "[";
	 	if(openChar.equals(")"))
	 		return "(";
	 	if(openChar.equals(">"))
	 		return "<";
	 	if(openChar.equals("}"))
	 		return "{";
	 	if(openChar.equals(":"))
	 		return ":";
	 	if(openChar.equals(";"))
	 		return ";";
	 	else
	 		return "";
	}

}
