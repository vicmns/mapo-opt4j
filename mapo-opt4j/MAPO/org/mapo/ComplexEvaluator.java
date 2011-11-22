package org.mapo;

import static org.opt4j.core.Objective.Sign.MIN;
import static org.opt4j.core.Objective.Sign.MAX;
import static org.opt4j.core.Objective.INFEASIBLE;

import java.util.ArrayList;

import org.opt4j.core.Objective;
import org.opt4j.core.Objectives;
import org.opt4j.core.problem.Evaluator;

import com.google.inject.Inject;

public class ComplexEvaluator implements Evaluator<ComplexPhenotype> {
	
	private double R=1.9872; //Gas Constant in cal/K-mol
	private double Ct=0.0001; //Molar Strand concentration, defualt 1M
	private double naClConct = 4.0; //Salt concentration
	private double nRatio = (4.35*Math.log10(naClConct)+3.5)*Math.pow(10, 5);
	private String mGene;
	private String mRNA;
	private ArrayList<Object> complex;
	private int mGeneLen;
	//private double dG=0.0, dH=0.0, dS=0.0;
	Objective DNAFreeEnergy = new Objective("DNAFreeEnergy", MIN);
	Objective RNAFreeEnergy = new Objective("RNAFreeEnergy", MIN);
	Objective DNAKf = new Objective("DNAKf", MAX);
	Objective DNAKr = new Objective("DNAKr", MIN);
	Objective RNAKr = new Objective("RNAKr", MIN);
	Objective length = new Objective("Length", MIN);
	
	@Inject
	public ComplexEvaluator(ComplexProblem problem){
		this.mGene=problem.getAB();
		this.mGeneLen=problem.getMGeneLen();
		this.mRNA=problem.getMRNA();
	}
		
	public Objectives evaluate(ComplexPhenotype phenotype) {
		this.complex=phenotype;
		Objectives obj = new Objectives();
		final int len= Integer.parseInt(complex.get(1).toString())-Integer.parseInt(complex.get(0).toString())+1;
		if(Integer.parseInt(complex.get(1).toString())>mGeneLen){
			obj.add(DNAFreeEnergy, INFEASIBLE);
			obj.setFeasible(false);
			obj.add(DNAKf,INFEASIBLE);
			obj.add(DNAKr,INFEASIBLE);
			obj.add(RNAKr,INFEASIBLE);
			obj.add(RNAFreeEnergy,INFEASIBLE);
			obj.add(length, INFEASIBLE);
		}
		else{
			double[] dnaEnergy=calculateEnergy(mGene.substring(Integer.parseInt(complex.get(0).toString()),
					Integer.parseInt(complex.get(1).toString())));
			double[] rnaEnergy=calculateRNADNAEnergy(this.mRNA.substring(Integer.parseInt(complex.get(0).toString()),
					Integer.parseInt(complex.get(1).toString())),true);
			if(dnaEnergy[3]<37.0){
				obj.add(DNAFreeEnergy, INFEASIBLE);
				obj.add(DNAKf,INFEASIBLE);
				obj.add(DNAKr,INFEASIBLE);
				obj.add(RNAKr,INFEASIBLE);
				obj.add(RNAFreeEnergy, INFEASIBLE);
				obj.add(length, INFEASIBLE);
				obj.setFeasible(false);
			}
			else{				
				obj.add(DNAFreeEnergy, dnaEnergy[2]);
				obj.add(DNAKf,dnaEnergy[4]);
				obj.add(DNAKr,dnaEnergy[5]);
				obj.add(RNAKr,rnaEnergy[5]);
				obj.add(RNAFreeEnergy, rnaEnergy[2]);
				obj.add(length, len);
			}
		}
		return obj;
	}
	
	private double[] cDuplexEnergy(String complex, boolean isInitial, boolean isTerminal){
		double dH=0,dS=0,dG=0;
		if(isInitial){
			dH+=0.2; dS+=-5.7; dG+=1.96; 
		}
		for(int i=0; i<complex.length()-1;i++){
			if((complex.charAt(i)=='A' && complex.charAt(i+1)=='A') || 
					(complex.charAt(i)=='T' && complex.charAt(i+1)=='T')){
				dH+=-7.6; dS+=-21.3; dG+=-1.0; 
			}
			if((complex.charAt(i)=='A' && complex.charAt(i+1)=='T')){
				dH+=-7.2; dS+=-20.4; dG+=-0.88; 
			}
			if((complex.charAt(i)=='T' && complex.charAt(i+1)=='A')){
				dH+=-7.2; dS+=-21.3; dG+=-0.58; 
			}
			if((complex.charAt(i)=='C' && complex.charAt(i+1)=='A') || 
					(complex.charAt(i)=='T' && complex.charAt(i+1)=='G')){
				dH+=-8.5; dS+=-22.7; dG+=-1.45;
			}
			if((complex.charAt(i)=='G' && complex.charAt(i+1)=='T')|| 
					(complex.charAt(i)=='A' && complex.charAt(i+1)=='C')){
				dH+=-8.4; dS+=-22.4; dG+=-1.44;
			}
			if((complex.charAt(i)=='C' && complex.charAt(i+1)=='T')| 
					(complex.charAt(i)=='A' && complex.charAt(i+1)=='G')){
				dH+=-7.8; dS+=-21.0; dG+=-1.28; 
			}
			if((complex.charAt(i)=='G' && complex.charAt(i+1)=='A') | 
					(complex.charAt(i)=='T' && complex.charAt(i+1)=='C')){
				dH+=-8.2; dS+=-22.2; dG+=-1.30; 
			}
			if((complex.charAt(i)=='C' && complex.charAt(i+1)=='G')){
				dH+=-10.6; dS+=-27.2; dG+=-2.17; 
			}
			if((complex.charAt(i)=='G' && complex.charAt(i+1)=='C')){
				dH+=-9.8; dS+=-24.4; dG+=-2.24; 
			}
			if((complex.charAt(i)=='G' && complex.charAt(i+1)=='G') || 
					(complex.charAt(i)=='C' && complex.charAt(i+1)=='C')){
				dH+=-8.0; dS+=-19.9; dG+=-1.84; 
			}
		}
		if(isTerminal && complex.charAt(complex.length()-1)=='A'){
			dG+=0.05;
		}
		//TODO: Add Symmetry correction?
		double[] pEnergy= {dH,dS,dG};//predicted Energy
		return pEnergy;
	}
	
	private double internalLoopGEnergy(int loopSize){
		double dG=0;
		if(loopSize==1)
			//Loop length not allowed
			//TODO internal mismatch
		if(loopSize==2)
			//TODO: Mismatch NN Parameters function
		if(loopSize==3)
			dG=3.2;
		if(loopSize==4)
			dG=3.6;
		if(loopSize==5)
			dG=4.0;
		if(loopSize==6)
			dG=4.4;
		if(loopSize==7)
			dG=4.6;
		if(loopSize==8)
			dG=4.8;
		if(loopSize==9)
			dG=4.9;
		if(loopSize==10 || loopSize==11)
			dG=4.9;
		if(loopSize==12 || loopSize==13)
			dG=5.2;
		if(loopSize==14 || loopSize==15)
			dG=5.4;
		if(loopSize==16 || loopSize==17)
			dG=5.6;
		if(loopSize==18 || loopSize==19)
			dG=5.8;
		if(loopSize==20 || loopSize<25)
			dG=5.9;
		if(loopSize==25 || loopSize<30)
			dG=6.3;
		if(loopSize==30)
			dG=6.6;
		if(loopSize>30){
			dG=6.6+2.44*R*310.15*Math.log(loopSize/30);
		}
		return dG;
	}
	
	private double[] fiveOverhangEnergy(char fiveEnd, char fiveFirstHang){
		double dH=0, dG=0;
		if(fiveEnd=='A'){
			if(fiveFirstHang=='A'){
				dH+=0.2; dG+=-0.51;
			}
			if(fiveFirstHang=='C'){
				dH+=0.6; dG+=-0.96;
			}
			if(fiveFirstHang=='G'){
				dH+=-1.1; dG+=-0.62;
			}
			if(fiveFirstHang=='T'){
				dH+=-6.9; dG+=-0.71;
			}
		}
		if(fiveEnd=='C'){
			if(fiveFirstHang=='A'){
				dH+=-6.3; dG+=-0.96;
			}
			if(fiveFirstHang=='C'){
				dH+=-4.4; dG+=-0.52;
			}
			if(fiveFirstHang=='G'){
				dH+=-5.1; dG+=-0.72;
			}
			if(fiveFirstHang=='T'){
				dH+=-4.0; dG+=-0.58;
			}
		}
		if(fiveEnd=='G'){
			if(fiveFirstHang=='A'){
				dH+=-3.7; dG+=-0.58;
			}
			if(fiveFirstHang=='C'){
				dH+=-4.0; dG+=-0.34;
			}
			if(fiveFirstHang=='G'){
				dH+=-3.9; dG+=-0.56;
			}
			if(fiveFirstHang=='T'){
				dH+=-4.9; dG+=-0.61;
			}
		}
		if(fiveEnd=='T'){
			if(fiveFirstHang=='A'){
				dH+=-2.9; dG+=-0.50;
			}
			if(fiveFirstHang=='C'){
				dH+=-4.1; dG+=-0.02;
			}
			if(fiveFirstHang=='G'){
				dH+=-4.2; dG+=0.48;
			}
			if(fiveFirstHang=='T'){
				dH+=-0.2; dG+=-0.10;
			}
		}
		double[] pEnergy={dH,dG};
		return pEnergy;
	}
	
	private double[] threeOverhangEnergy(char threeEnd, char threeFirstHang){
		double dH=0, dG=0;
		//Three end
		if(threeEnd=='A'){
			if(threeFirstHang=='A'){
				dH+=-0.5; dG+=-0.12;
			}
			if(threeFirstHang=='C'){
				dH+=4.7; dG+=0.28;
			}
			if(threeFirstHang=='G'){
				dH+=-4.1; dG+=-0.01;
			}
			if(threeFirstHang=='T'){
				dH+=-3.8; dG+=0.13;
			}
		}
		if(threeEnd=='C'){
			if(threeFirstHang=='A'){
				dH+=-5.9; dG+=-0.82;
			}
			if(threeFirstHang=='C'){
				dH+=-2.6; dG+=-0.31;
			}
			if(threeFirstHang=='G'){
				dH+=-3.2; dG+=-0.01;
			}
			if(threeFirstHang=='T'){
				dH+=-5.2; dG+=-0.52;
			}
		}
		if(threeEnd=='G'){
			if(threeFirstHang=='A'){
				dH+=-2.1; dG+=-0.92;
			}
			if(threeFirstHang=='C'){
				dH+=4.7; dG+=0.28;
			}
			if(threeFirstHang=='G'){
				dH+=-4.1; dG+=-0.01;
			}
			if(threeFirstHang=='T'){
				dH+=-3.8; dG+=0.13;
			}
		}
		if(threeEnd=='T'){
			if(threeFirstHang=='A'){
				dH+=-0.7; dG+=-0.48;
			}
			if(threeFirstHang=='C'){
				dH+=4.4; dG+=-0.19;
			}
			if(threeFirstHang=='G'){
				dH+=-1.6; dG+=-0.50;
			}
			if(threeFirstHang=='T'){
				dH+=2.9; dG+=-0.29;
			}
		}
		double[] pEnergy={dH,dG};
		return pEnergy;
	}
	
	private double bulgesEnergy(int loopSize){
		double dG=0;
		if(loopSize==1)
			dG=4.0;
		//TODO: intervening base pair stack
		if(loopSize==2)
			dG=2.9;
		if(loopSize==3)
			dG=3.1;
		if(loopSize==4)
			dG=3.2;
		if(loopSize==5)
			dG=3.3;
		if(loopSize==6)
			dG=3.5;
		if(loopSize==7)
			dG=3.7;
		if(loopSize==8)
			dG=3.9;
		if(loopSize==9)
			dG=4.1;
		if(loopSize==10 || loopSize==11)
			dG=4.3;
		if(loopSize==12 || loopSize==13)
			dG=4.5;
		if(loopSize==14 || loopSize==15)
			dG=4.6;
		if(loopSize==16 || loopSize==17)
			dG=5.0;
		if(loopSize==18 || loopSize==19)
			dG=5.2;
		if(loopSize==20 || loopSize<25)
			dG=5.3;
		if(loopSize==25 || loopSize<30)
			dG=5.6;
		if(loopSize==30)
			dG=5.9;
		if(loopSize>30)
			dG=5.9+2.44*R*310.15*Math.log(loopSize/30);
		return dG;
	}
	
	public double[] calculateEnergy(String mGene){
		double dH=0,dS=0,dG=0;
		int idx=0;
		boolean isInitial=true;
		boolean isTerminal=false;
		int numPairBases=0;
		int totalOverhang=0;
		for(int i=3; i<complex.size(); i+=3){
			//Overhang
			if(complex.get(i).equals(":")){
				idx+=Integer.parseInt(complex.get(i+1).toString());
				totalOverhang+=Integer.parseInt(complex.get(i+1).toString());
				double[] energy = fiveOverhangEnergy(mGene.charAt(idx-1), mGene.charAt(idx));
				dH+=energy[0]; dG+=energy[1];
				isInitial=false;
			}
			if(complex.get(i).equals("[")){
				if(i>=complex.size()-1){
					isTerminal=true;
				}
				double[] energy = cDuplexEnergy(mGene.substring(idx, idx+Integer.parseInt(complex.get(i+1).toString()) ), isInitial, isTerminal);
				idx+=Integer.parseInt(complex.get(i+1).toString());
				numPairBases+=Integer.parseInt(complex.get(i+1).toString());
				dH+=energy[0]; dS+=energy[1]; dG+=energy[2];
				isInitial=false;
			}
			if(complex.get(i).equals("(")){
				double energy=internalLoopGEnergy(Integer.parseInt(complex.get(i+1).toString()));
				dG+=energy;
				idx+=Integer.parseInt(complex.get(i+1).toString());
				isInitial=false;
			}
			if(complex.get(i).equals("{") || complex.get(i).equals("<")){
				double energy=bulgesEnergy(Integer.parseInt(complex.get(i+1).toString()));
				dG+=energy;
				idx+=Integer.parseInt(complex.get(i+1).toString());
				totalOverhang+=Integer.parseInt(complex.get(i+1).toString());
				isInitial=false;
			}
			if(complex.get(i).equals(";")){
				double[] energy = threeOverhangEnergy(mGene.charAt(idx-1), mGene.charAt(idx));
				idx+=Integer.parseInt(complex.get(i+1).toString());
				dH+=energy[0]; dG+=energy[1];
			}
		}
		double Tm=dH*1000/(dS+R*Math.log(Ct/4))-273.15;
		//Calculate Kinetics
		int shortestLen = Integer.parseInt(complex.get(1).toString())-Integer.parseInt(
				complex.get(0).toString())-totalOverhang+1;
		double[] kinetics = calculateKinetics(shortestLen, numPairBases, dG, Tm);
		double[] pEnergy= {dH,dS,dG,Tm,kinetics[0],kinetics[1]};//predicted Energy
		return pEnergy;
	}
	
	/**
	 * Calculate the thermodynamics parameters of RNA/DNA complexes
	 * @param complex
	 * 			The RNA complex to calculate thermodynamics parameters
	 * 			Note: The parameters are calculated implying complete base-pairing  .
	 * @param isInitial
	 * 			Boolean variable that defines if the complex to calculate is initial or not.
	 */
	public double[] calculateRNADNAEnergy(String mRNA, boolean isInitial){
		double dH=0,dS=0,dG=0;
		if(isInitial){
			dH+=1.9; dS+=-3.9; dG+=3.1; 
		}
		for(int i=0; i<mRNA.length()-1;i++){
			if((mRNA.charAt(i)=='A' && mRNA.charAt(i+1)=='A')){
				dH+=-7.8; dS+=-21.9; dG+=-1.0; 
			}
			if((mRNA.charAt(i)=='A' && mRNA.charAt(i+1)=='C')){
				dH+=-5.9; dS+=-12.3; dG+=-2.1; 
			}
			if((mRNA.charAt(i)=='A' && mRNA.charAt(i+1)=='G')){
				dH+=-9.1; dS+=-23.5; dG+=-1.8; 
			}
			if((mRNA.charAt(i)=='A' && mRNA.charAt(i+1)=='U')){
				dH+=-8.3; dS+=-23.9; dG+=-0.9;
			}
			if((mRNA.charAt(i)=='C' && mRNA.charAt(i+1)=='A')){
				dH+=-9.0; dS+=-26.1; dG+=-0.9;
			}
			if((mRNA.charAt(i)=='C' && mRNA.charAt(i+1)=='C')){
				dH+=-9.3; dS+=-23.2; dG+=-2.1; 
			}
			if((mRNA.charAt(i)=='C' && mRNA.charAt(i+1)=='G')){
				dH+=-16.3; dS+=-47.1; dG+=-1.7; 
			}
			if((mRNA.charAt(i)=='C' && mRNA.charAt(i+1)=='U')){
				dH+=-7.0; dS+=-19.7; dG+=-0.9; 
			}
			if((mRNA.charAt(i)=='G' && mRNA.charAt(i+1)=='A')){
				dH+=-5.5; dS+=-13.5; dG+=-1.3; 
			}
			if((mRNA.charAt(i)=='G' && mRNA.charAt(i+1)=='C')){
				dH+=-8.0; dS+=-17.1; dG+=-2.7; 
			}
			if((mRNA.charAt(i)=='G' && mRNA.charAt(i+1)=='G')){
				dH+=-12.8; dS+=-31.9; dG+=-2.9; 
			}
			if((mRNA.charAt(i)=='G' && mRNA.charAt(i+1)=='U')){
				dH+=-7.8; dS+=-21.6; dG+=-1.1; 
			}
			if((mRNA.charAt(i)=='U' && mRNA.charAt(i+1)=='A')){
				dH+=-7.8; dS+=-23.2; dG+=-0.6; 
			}
			if((mRNA.charAt(i)=='U' && mRNA.charAt(i+1)=='C')){
				dH+=-8.6; dS+=-22.9; dG+=-1.5; 
			}
			if((mRNA.charAt(i)=='U' && mRNA.charAt(i+1)=='G')){
				dH+=-10.4; dS+=-28.4; dG+=-1.6; 
			}
			if((mRNA.charAt(i)=='U' && mRNA.charAt(i+1)=='U')){
				dH+=-11.5; dS+=-36.4; dG+=-0.2; 
			}
		}
		double Tm=dH*1000/(dS+R*Math.log(Ct/4))-273.15;
		double[] kinetics = calculateKinetics(mRNA.length()-1, mRNA.length()-1, dG, Tm);
		double[] pEnergy= {dH,dS,dG,Tm,kinetics[0],kinetics[1]};//predicted Energy
		return pEnergy;
	}
	
	private double[] calculateKinetics(int shortestLen, int numPairBases, double dG, double Tm){
		double incTemp=(Tm+273.15)-298.15;
		double kF=(this.nRatio*Math.sqrt(shortestLen))/numPairBases;
		double kR=kF*Math.exp(dG/(this.R*incTemp));
		double[] kinectics = {kF, kR};
		return kinectics;
	}
}
