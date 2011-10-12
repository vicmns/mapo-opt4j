package org.mapo;

import java.util.ArrayList;
import java.util.Random;

import org.opt4j.common.random.Rand;
import org.opt4j.genotype.ListGenotype;
import org.opt4j.genotype.StringGenotype;
import org.opt4j.operator.crossover.Crossover;
import org.opt4j.optimizer.ea.Pair;

public class CrossOverComplex<G extends StringGenotype<?>> implements Crossover<G> {
	private Random random;
	private Object[] p1SizeParse;
	private Object[] p1SignParse;
	private Object[] p2SizeParse;
	private Object[] p2SignParse;
	private int idxSp1;
	private int idxEp1;
	private int alphaP1;
	private int idxSp2;
	private int idxEp2;
	private int alphaP2;
	
	
	public CrossOverComplex(Rand random) {
		this.random = random;
	}
	
	@Override
	public  Pair<G> crossover(G p1, G p2) {
		Object[][] p1Parse=complexParse(p1.toString());
		this.p1SizeParse=p1Parse[0];
		this.p1SignParse=p1Parse[1];
		this.idxSp1=Integer.parseInt(p1Parse[2][0].toString());
		this.idxEp1=Integer.parseInt(p1Parse[2][1].toString());
		this.alphaP1=Integer.parseInt(p1Parse[2][2].toString());
		Object[][] p2Parse=complexParse(p2.toString());
		this.p2SizeParse=p2Parse[0];
		this.p2SignParse=p2Parse[1];
		this.idxSp2=Integer.parseInt(p2Parse[2][0].toString());
		this.idxEp2=Integer.parseInt(p2Parse[2][1].toString());
		this.alphaP2=Integer.parseInt(p2Parse[2][2].toString());
		//New Instances for childs
		ListGenotype<Object> o1 = p1.newInstance();
		ListGenotype<Object> o2 = p2.newInstance();
		//HERE
		ArrayList<String> child1 = new ArrayList<String>();
		ArrayList<String> child2 = new ArrayList<String>();
		//String[] childs = new String[2];
		int fAlphaP1=this.idxSp1;
		int fAlphaP2=this.idxSp2;
		int idxMP1=0; //Mutation index on parent 1
		int idxMP2=0; //Mutation index on parent 2
		int idxSC1=this.idxSp1;
		int idxSC2=this.idxSp2;
		int idxEC1=this.idxEp1;
		int idxEC2=this.idxEp2;
		int alphaC1=this.alphaP1;
		int alphaC2=this.alphaP2;
		int i=0;
		/*
		 * Find where the mutation starts
		 */
		while(fAlphaP1!=this.alphaP1){
			fAlphaP1+=Integer.parseInt(p1SizeParse[i].toString());
			i++;
			idxMP1=i;
		}
		i=0;
		while(fAlphaP2!=this.alphaP2){
			fAlphaP2+=Integer.parseInt(p2SizeParse[i].toString());
			i++;
			idxMP2=i;
		}
		if(random.nextBoolean()){
			int sumIdx=0;
			//Take the left configuration before mutation from parent 1 and left configuration (after) from parent 2
			for(int j=0; j<idxMP1; j++){
				sumIdx+=Integer.parseInt(p1SizeParse[j].toString());
				child1.add(getOpenChar(p1SignParse[j].toString())+p1SizeParse[j].toString()+p1SignParse[j].toString());
			}
			idxEC1=sumIdx;
			sumIdx=0;
			for(int j=0; j<idxMP2;j++){
				sumIdx+=Integer.parseInt(p2SizeParse[j].toString());
				child2.add(getOpenChar(p2SignParse[j].toString())+p2SizeParse[j].toString()+p2SignParse[j].toString());
			}
			idxEC2=sumIdx;
			sumIdx=0;
			child1.add(getOpenChar(p1SignParse[idxMP1].toString())+p1SizeParse[idxMP1].toString()+p1SignParse[idxMP1].toString());
			idxEC1+=Integer.parseInt(p1SizeParse[idxMP1].toString());
			child2.add(getOpenChar(p2SignParse[idxMP2].toString())+p2SizeParse[idxMP2].toString()+p2SignParse[idxMP2].toString());
			idxEC2+=Integer.parseInt(p2SizeParse[idxMP2].toString());
			for(int j=idxMP2+1; j<p2SignParse.length; j++){
				sumIdx+=Integer.parseInt(p2SizeParse[j].toString());
				child1.add(getOpenChar(p2SignParse[j].toString())+p2SizeParse[j].toString()+p2SignParse[j].toString());
			}
			idxEC1+=(idxSC1+sumIdx);
			sumIdx=0;
			for(int j=idxMP1+1; j<p1SignParse.length; j++){
				sumIdx+=Integer.parseInt(p1SizeParse[j].toString());
				child2.add(getOpenChar(p1SignParse[j].toString())+p1SizeParse[j].toString()+p1SignParse[j].toString());
			}
			idxEC2+=(idxSC2+sumIdx);
			child1.add(0,idxSC1+","+idxEC1+","+alphaC1);
			child2.add(0,idxSC2+","+idxEC2+","+alphaC2);
		}
		else{
			idxSC1=this.idxSp2;
			idxSC2=this.idxSp1;
			int sumAlpha=0;
			int sumIdx=0;
			//Take the left configuration from parent 2 and left configuration from parent 1
			//Recalculate Ending indexes and alpha
			for(int j=0; j<idxMP2;j++){
				sumIdx+=Integer.parseInt(p2SizeParse[j].toString());
				sumAlpha+=Integer.parseInt(p2SizeParse[j].toString());
				child1.add(getOpenChar(p2SignParse[j].toString())+p2SizeParse[j].toString()+p2SignParse[j].toString());
			}
			idxEC1=sumIdx;
			alphaC1=sumAlpha;
			sumAlpha=0;
			sumIdx=0;
			for(int j=0; j<idxMP1; j++){
				sumIdx+=Integer.parseInt(p1SizeParse[j].toString());
				sumAlpha+=Integer.parseInt(p1SizeParse[j].toString());
				child2.add(getOpenChar(p1SignParse[j].toString())+p1SizeParse[j].toString()+p1SignParse[j].toString());
			}
			idxEC2=sumIdx;
			alphaC2=sumAlpha;
			sumIdx=0;
			child1.add(getOpenChar(p1SignParse[idxMP1].toString())+p1SizeParse[idxMP1].toString()+p1SignParse[idxMP1].toString());
			idxEC1+=Integer.parseInt(p1SizeParse[idxMP1].toString());
			child2.add(getOpenChar(p2SignParse[idxMP2].toString())+p2SizeParse[idxMP2].toString()+p2SignParse[idxMP2].toString());
			idxEC2+=Integer.parseInt(p2SizeParse[idxMP2].toString());
			for(int j=idxMP1+1; j<p1SignParse.length; j++){
				sumIdx+=Integer.parseInt(p1SizeParse[j].toString());
				child1.add(getOpenChar(p1SignParse[j].toString())+p1SizeParse[j].toString()+p1SignParse[j].toString());
			}
			idxEC1+=(idxSC1+sumIdx);
			sumIdx=0;
			for(int j=idxMP2+1; j<p2SignParse.length; j++){
				sumIdx+=Integer.parseInt(p2SizeParse[j].toString());
				child2.add(getOpenChar(p2SignParse[j].toString())+p2SizeParse[j].toString()+p2SignParse[j].toString());
			}
			idxEC2+=(idxSC2+sumIdx);
			child1.add(0,idxSC1+","+idxEC1+","+(alphaC1+idxSC1));
			child2.add(0,idxSC2+","+idxEC2+","+(alphaC2+idxSC2));
		}
		//childs[0]=toString(child1.toArray());
		//childs[1]=toString(child2.toArray());
		o1.addAll(child1);
		o2.addAll(child2);
		//HERE
		@SuppressWarnings("unchecked")
		Pair<G> offspring = new Pair<G>((G) o1, (G) o2);
		return offspring;
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
	
	@SuppressWarnings("unused")
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
