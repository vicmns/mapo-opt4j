package org.mapo;


import com.google.inject.Inject;

public class ComplexProblem {

	protected int idxS;
	protected int idxE;
	protected String Gene; //Correct cDNA/DNA Sequence
	protected String AB; //Mutated cDNA/RNA
	
	@Inject
	public ComplexProblem(){
		String gene, ab;
		gene=readFromFile();
		if(gene=="")
			return;
		ab=readFromFile();
		if(ab=="")
			return;
		this.Gene=gene;
		this.AB=ab;
		int[] idx=findMutation(this.Gene, this.AB);
		this.idxS=idx[0];
		this.idxE=idx[1];
	}
	
	public int getIdxS(){
		return this.idxS;
	}
	
	public int getIdxE(){
		return this.idxE;
	}
	
	public int getGeneLength(){
		return Gene.length();
	}
	
	public int getABLength(){
		return AB.length();
	}
	
	public String getGene(){
		return Gene;
	}
	
	public String getAB(){
		return AB;
	}
	
	private String readFromFile(){
		FileManager fileMan= new FileManager();
		return fileMan.readStringFromFile();
	}
	
	private  int[] findMutation(String a, String b){
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
}
