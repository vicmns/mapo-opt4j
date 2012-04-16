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
	
	private double R=1.987; //Gas Constant in cal/K-mol
	private double Ct=0.000001; //Molar Strand concentration, defualt 1uM
	private double naClConct = 4.0; //Salt concentration
	private double nRatio = (4.35*Math.log10(naClConct)+3.5)*Math.pow(10, 5);
	private String mGene;
	private String mRNA;
	private String Gene;
	private String h_mRNA;
	private ArrayList<Object> complex;
	private int mGeneLen;
	private int numPairBases;
	private int totalOverHang;
	private int len;
	//private double dG=0.0, dH=0.0, dS=0.0;
	Objective DNAFreeEnergy = new Objective("DNAFreeEnergy", MIN);
	Objective RNAFreeEnergy = new Objective("RNAFreeEnergy", MIN);
	Objective hRNAFreeEnergy = new Objective("hRNAFreeEnergy", MAX);
	Objective DNAKf = new Objective("DNAKf", MAX);
	Objective DNAKr = new Objective("DNAKr", MIN);
	Objective RNAKf = new Objective("RNAKf", MAX);
	Objective RNAKr = new Objective("RNAKr", MIN);
	Objective hRNAKf = new Objective("hRNAKf", MIN);
	Objective hRNAKr = new Objective("hRNAKr", MAX);
	Objective length = new Objective("Length", MIN);
	
	@Inject
	public ComplexEvaluator(ComplexProblem problem){
		this.mGene=problem.getAB();
		this.mGeneLen=problem.getMGeneLen();
		this.mRNA=problem.getMRNA();
		this.Gene=problem.getGene();
		this.h_mRNA=problem.getHmRNA();
	}
		
	public Objectives evaluate(ComplexPhenotype phenotype) {
		int idxInd = 0;
		int i=4;
		int sumNuc=0;
		this.complex=phenotype;
		while(Integer.parseInt(phenotype.get(2).toString())!=sumNuc){
			if(!phenotype.get(i-1).equals("{"))
				sumNuc+=Integer.parseInt(phenotype.get(i).toString());
			i+=3;
		}
		idxInd=i+1;
		if(!phenotype.get(idxInd).equals(";") && !phenotype.get(idxInd).equals(":") && !phenotype.get(idxInd).equals(")")){
			System.out.println("HAHA! theres a problem!");
			System.out.println(phenotype.toString());
		}
		Objectives obj = new Objectives();
		this.len = Integer.parseInt(complex.get(1).toString()) - Integer.parseInt(complex.get(0).toString()) + 1;
		if(Integer.parseInt(complex.get(1).toString()) + 1 > mGeneLen || Integer.parseInt(complex.get(0).toString()) < 0){
			obj.add(DNAFreeEnergy, INFEASIBLE);
			obj.setFeasible(false);
			obj.add(DNAKf,INFEASIBLE);
			obj.add(DNAKr,INFEASIBLE);
			obj.add(RNAFreeEnergy,INFEASIBLE);
			obj.add(RNAKf,INFEASIBLE);
			obj.add(RNAKr,INFEASIBLE);
			obj.add(hRNAFreeEnergy,INFEASIBLE);
			obj.add(hRNAKf,INFEASIBLE);
			obj.add(hRNAKr,INFEASIBLE);
			obj.add(length, INFEASIBLE);
		}
		else{
	        boolean isInitial = true;
			double[] dnaEnergy=calculateDNAEnergy(this.mGene.substring(Integer.parseInt(complex.get(0).toString()),
					Integer.parseInt(complex.get(1).toString()) + 1));
			/*double[] rna_dnaEnergy= calculateRNADNAKinectics(this.mRNA.substring(Integer.parseInt(complex.get(0).toString()),
					Integer.parseInt(complex.get(1).toString())),true);*/
			double[] rnaEnergy=calculateRNADNAEnergy(this.mRNA.substring(Integer.parseInt(complex.get(0).toString()),
					Integer.parseInt(complex.get(1).toString()) + 1),true);
			double[] hRnaEnergy = calculateHealthyRna_DxEnergy();
			if(dnaEnergy[3]<37.0 || dnaEnergy[4] / dnaEnergy[5] < 1.5){
				obj.add(DNAFreeEnergy, INFEASIBLE);
				obj.add(DNAKf,INFEASIBLE);
				obj.add(DNAKr,INFEASIBLE);
				obj.add(RNAFreeEnergy, INFEASIBLE);
				obj.add(RNAKf,INFEASIBLE);
				obj.add(RNAKr,INFEASIBLE);
				obj.add(hRNAFreeEnergy,INFEASIBLE);
				obj.add(hRNAKf,INFEASIBLE);
				obj.add(hRNAKr,INFEASIBLE);
				obj.add(length, INFEASIBLE);
				obj.setFeasible(false);
			}
			else{				
				obj.add(DNAFreeEnergy, dnaEnergy[2]);
				obj.add(DNAKf,dnaEnergy[4]);
				obj.add(DNAKr,dnaEnergy[5]);
				obj.add(RNAFreeEnergy, rnaEnergy[2]);
				obj.add(RNAKf,rnaEnergy[4]);
				obj.add(RNAKr,rnaEnergy[5]);
				obj.add(hRNAFreeEnergy,hRnaEnergy[2]);
				obj.add(hRNAKf,hRnaEnergy[5]);
				obj.add(hRNAKr,hRnaEnergy[4]);
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
		if(isTerminal && (complex.charAt(complex.length()-1)=='A' || 
				complex.charAt(complex.length()-1)=='T')){
			dH+=2.2; dS+=6.9; dG+=0.05;
		}
		double[] pEnergy= {dH,dS,dG};//predicted Energy
		return pEnergy;
	}
	
	private double[] cDuplexRNAEnergy(String mRNA, boolean isInitial, boolean isTerminal){
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
		double energy[] = {dH,dS,dG};
		return energy;
	}
	
	private double internalLoopGEnergy(int loopSize){
		double dG=0;
		if(loopSize==2)
			dG=3.2;
		else if(loopSize==3)
			dG=3.2;
		else if(loopSize==4)
			dG=3.6;
		else if(loopSize==5)
			dG=4.0;
		else if(loopSize==6)
			dG=4.4;
		else if(loopSize==7)
			dG=4.6;
		else if(loopSize==8)
			dG=4.8;
		else if(loopSize==9)
			dG=4.9;
		else if(loopSize==10 || loopSize==11)
			dG=4.9;
		else if(loopSize==12 || loopSize==13)
			dG=5.2;
		else if(loopSize==14 || loopSize==15)
			dG=5.4;
		else if(loopSize==16 || loopSize==17)
			dG=5.6;
		else if(loopSize==18 || loopSize==19)
			dG=5.8;
		else if(loopSize==20 || loopSize<25)
			dG=5.9;
		else if(loopSize==25 || loopSize<30)
			dG=6.3;
		else if(loopSize==30)
			dG=6.6;
		else if(loopSize>30){
			dG=6.6+2.44*R*310.15*Math.log(loopSize/30);
		}
		return dG;
	}
	
	private double internalRNALoopGEnergy(int loopSize) {
        double dG = 0;
        if (loopSize == 2) { 
            dG=0.8;
        }
        else if (loopSize == 3) {
            dG = 1.3;
        }
        else if (loopSize == 4) {
            dG = 1.7;
        }
        else if (loopSize == 5) {
            dG = 2.1;
        }
        else if (loopSize == 6) {
            dG = 2.5;
        }
        else if (loopSize == 7) {
            dG = 2.6;
        }
        else if (loopSize == 8) {
            dG = 2.8;
        }
        else if (loopSize == 9) {
            dG = 3.1;
        }
        else if (loopSize == 10 || loopSize == 11) {
            dG = 3.6;
        }
        else if (loopSize == 12 || loopSize == 13) {
            dG = 4.4;
        }
        else if (loopSize == 14 || loopSize == 15) {
            dG = 5.1;
        }
        else if (loopSize == 16 || loopSize == 17) {
            dG = 5.6;
        }
        else if (loopSize == 18 || loopSize == 19) {
            dG = 6.2;
        }
        else if (loopSize == 20 || loopSize < 25) {
            dG = 6.6;
        }
        else if (loopSize == 25 || loopSize < 30) {
            dG = 7.6;
        }
        else if (loopSize == 30) {
            dG = 8.4;
        }
        else if (loopSize > 30) {
            dG = 8.4 + 2.44 * R * 310.15 * Math.log(loopSize / 30);
        }
        return dG;
    }
	
	private double internalSingleMissRNA(String mGene) {
        double energy = 0;
        //Trimmers energy
        int j=1;
        for (int i = 0; i < mGene.length()-1; i += 2) {
            if(i==2){
                j++;
                i--;
            }
            if (mGene.charAt(j) == 'G') {
                //AA configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'A') {
                    energy += 0.81;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'A') {
                    energy += 0.79;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'U') {
                    energy += 1.28;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'U') {
                    energy += 1.41;
                }
                //GC configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += -0.91;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += -0.19;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += -0.68;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += -1.40;
                }
                //AG configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'G') {
                    energy += -0.02;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'G') {
                    energy += 0.11;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'C') {
                    energy += -0.38;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'G') {
                    energy += -0.02;
                }
                //GA configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += -0.08;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += 0.64;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += 0.39;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += 1.11;
                }
            } else if (mGene.charAt(j) == 'A') {
                //AA configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'A') {
                    energy += 2.43;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'A') {
                    energy += 2.03;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'U') {
                    energy += 2.92;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'U') {
                    energy += 2.98;
                }
                //GC configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += 0.72;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += 1.11;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += 1.09;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += 0.70;
                }
                //AG configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'G') {
                    energy += 1.28;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'G') {
                    energy += 1.34;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'C') {
                    energy += 1.32;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'G') {
                    energy += 1.28;
                }
                //GA configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += 1.87;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += 2.26;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += 2.36;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += 2.75;
                }
            } else if (mGene.charAt(j) == 'U') {
                //AA configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'A') {
                    energy += 1.84;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'A') {
                    energy += 2.28;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'U') {
                    energy += 1.66;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'U') {
                    energy += 2.10;
                }
                //GG configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += 0.32;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += 0.63;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += 0.47;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += 0.16;
                }
                //AG configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'G') {
                    energy += 0.77;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'G') {
                    energy += 1.21;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'C') {
                    energy += 1.05;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'G') {
                    energy += 0.77;
                }
                //GA configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += 1.39;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += 1.70;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += 1.21;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += 1.52;
                }
            } else if (mGene.charAt(j) == 'C') {
                //AA configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'A') {
                    energy += 3.33;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'A') {
                    energy += 2.85;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'U') {
                    energy +=3.51;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'U') {
                    energy += 3.03;
                }
                //GC configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += 1.42;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += 1.50;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += 1.77;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += 1.69;
                }
                //AG configuration
                if (mGene.charAt(i) == 'A' && mGene.charAt(i+2) == 'G') {
                    energy += 2.09;
                }
                if (mGene.charAt(i) == 'U' && mGene.charAt(i+2) == 'G') {
                    energy += 1.61;
                }
                if (mGene.charAt(i) == 'U'&& mGene.charAt(i+2) == 'C') {
                    energy += 1.88;
                }
                if (mGene.charAt(i) == 'A'&& mGene.charAt(i+2) == 'G') {
                    energy += 2.09;
                }
                //GA configuration
                if (mGene.charAt(i) == 'G' && mGene.charAt(i+2) == 'G') {
                    energy += 2.66;
                }
                if (mGene.charAt(i) == 'C' && mGene.charAt(i+2) == 'G') {
                    energy += 2.74;
                }
                if (mGene.charAt(i) == 'G'&& mGene.charAt(i+2) == 'G') {
                    energy += 2.84;
                }
                if (mGene.charAt(i) == 'C'&& mGene.charAt(i+2) == 'C') {
                    energy += 2.92;
                }
            }
        }
        return energy;
    }
	
	private double[] fiveOverhangEnergy(char fiveFirstHang, char fiveEnd){
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
	
	private double[] fiveOverhangRNAEnergy(char fiveEnd, char fiveFirstHang) {
        double dH = 0, dG = 0;
        if (fiveEnd == 'A') {
            if (fiveFirstHang == 'A') {
                dH += 0.2;
                dG += -0.3;
            }
            if (fiveFirstHang == 'C') {
                dH += 0.6;
                dG += -0.3;
            }
            if (fiveFirstHang == 'G') {
                dH += -1.1;
                dG += -0.4;
            }
            if (fiveFirstHang == 'U') {
                dH += -6.9;
                dG += -0.2;
            }
        }
        if (fiveEnd == 'C') {
            if (fiveFirstHang == 'A') {
                dH += -6.3;
                dG += -0.5;
            }
            if (fiveFirstHang == 'C') {
                dH += -4.4;
                dG += -0.2;
            }
            if (fiveFirstHang == 'G') {
                dH += -5.1;
                dG += -0.2;
            }
            if (fiveFirstHang == 'U') {
                dH += -4.0;
                dG += -0.1;
            }
        }
        if (fiveEnd == 'G') {
            if (fiveFirstHang == 'A') {
                dH += -3.7;
                dG += -0.2;
            }
            if (fiveFirstHang == 'C') {
                dH += -4.0;
                dG += -0.3;
            }
            if (fiveFirstHang == 'G') {
                dH += -3.9;
                dG += -0.0;
            }
            if (fiveFirstHang == 'U') {
                dH += -4.9;
                dG += -0.0;
            }
        }
        if (fiveEnd == 'U') {
            if (fiveFirstHang == 'A') {
                dH += -2.9;
                dG += -0.3;
            }
            if (fiveFirstHang == 'C') {
                dH += -4.1;
                dG += -0.2;
            }
            if (fiveFirstHang == 'G') {
                dH += -4.2;
                dG += -0.2;
            }
            if (fiveFirstHang == 'U') {
                dH += -0.2;
                dG += -0.2;
            }
        }
        double[] pEnergy = {dH, dG};
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
	
	private double[] threeOverhangRNAEnergy(char threeEnd, char threeFirstHang){
        double dH = 0, dG = 0;
        //Three end
        if (threeEnd == 'A') {
            if (threeFirstHang == 'A') {
                dH += -0.5;
                dG += -0.8;
            }
            if (threeFirstHang == 'C') {
                dH += 4.7;
                dG += -0.5;
            }
            if (threeFirstHang == 'G') {
                dH += -4.1;
                dG += -0.8;
            }
            if (threeFirstHang == 'U') {
                dH += -3.8;
                dG += 0.6;
            }
        }
        if (threeEnd == 'C') {
            if (threeFirstHang == 'A') {
                dH += -5.9;
                dG += -1.7;
            }
            if (threeFirstHang == 'C') {
                dH += -2.6;
                dG += -0.8;
            }
            if (threeFirstHang == 'G') {
                dH += -3.2;
                dG += -1.7;
            }
            if (threeFirstHang == 'U') {
                dH += -5.2;
                dG += -1.2;
            }
        }
        if (threeEnd == 'G') {
            if (threeFirstHang == 'A') {
                dH += -2.1;
                dG += -1.1;
            }
            if (threeFirstHang == 'C') {
                dH += 4.7;
                dG += -0.4;
            }
            if (threeFirstHang == 'G') {
                dH += -4.1;
                dG += -1.3;
            }
            if (threeFirstHang == 'U') {
                dH += -3.8;
                dG += 0.6;
            }
        }
        if (threeEnd == 'U') {
            if (threeFirstHang == 'A') {
                dH += -0.7;
                dG += -0.7;
            }
            if (threeFirstHang == 'C') {
                dH += 4.4;
                dG += -0.1;
            }
            if (threeFirstHang == 'G') {
                dH += -1.6;
                dG += -0.7;
            }
            if (threeFirstHang == 'U') {
                dH += 2.9;
                dG += -0.1;
            }
        }
        double[] pEnergy = {dH, dG};
        return pEnergy;
    }
	
	private double bulgesEnergy(int loopSize){
		double dG=0;
		if(loopSize==1)
			dG=4.0;
		//TODO: intervening base pair stack
		if(loopSize==2)
			dG=2.9;
		else if(loopSize==3)
			dG=3.1;
		else if(loopSize==4)
			dG=3.2;
		else if(loopSize==5)
			dG=3.3;
		else if(loopSize==6)
			dG=3.5;
		else if(loopSize==7)
			dG=3.7;
		else if(loopSize==8)
			dG=3.9;
		else if(loopSize==9)
			dG=4.1;
		else if(loopSize==10 || loopSize==11)
			dG=4.3;
		else if(loopSize==12 || loopSize==13)
			dG=4.5;
		else if(loopSize==14 || loopSize==15)
			dG=4.6;
		else if(loopSize==16 || loopSize==17)
			dG=5.0;
		else if(loopSize==18 || loopSize==19)
			dG=5.2;
		else if(loopSize==20 || loopSize<25)
			dG=5.3;
		else if(loopSize==25 || loopSize<30)
			dG=5.6;
		else if(loopSize==30)
			dG=5.9;
		else if(loopSize>30)
			dG=5.9+2.44*R*310.15*Math.log(loopSize/30);
		return dG;
	}
	
	private double internalSingleMiss(String mGene){
		double energy = 0;
		int j=1;
        int i=0;
        for (int k=0; k<mGene.length(); k++ ) {
            if (mGene.charAt(i) == 'G') {
                if (mGene.charAt(j) == 'A') {
                    energy += 0.17;
                }
                if (mGene.charAt(j) == 'C') {
                    energy += 0.79;
                }
                if (mGene.charAt(j) == 'G') {
                    energy += -1.11;
                }
                if (mGene.charAt(j) == 'T') {
                    energy += 0.45;
                }
            } else if (mGene.charAt(i) == 'C') {
                if (mGene.charAt(j) == 'A') {
                    energy += 0.43;
                }
                if (mGene.charAt(j) == 'C') {
                    energy += 0.70;
                }
                if (mGene.charAt(j) == 'G') {
                    energy += -0.11;
                }
                if (mGene.charAt(j) == 'T') {
                    energy += -0.12;
                }
            } else if (mGene.charAt(i) == 'A') {
                if (mGene.charAt(j) == 'A') {
                    energy += 0.61;
                }
                if (mGene.charAt(j) == 'C') {
                    energy += 1.33;
                }
                if (mGene.charAt(j) == 'G') {
                    energy += -0.13;
                }
                if (mGene.charAt(j) == 'T') {
                    energy += 0.69;
                }
            } else if (mGene.charAt(i) == 'T') {
                if (mGene.charAt(j) == 'A') {
                    energy += 0.69;
                }
                if (mGene.charAt(j) == 'C') {
                    energy += 1.05;
                }
                if (mGene.charAt(j) == 'G') {
                    energy += 0.44;
                }
                if (mGene.charAt(j) == 'T') {
                    energy += 0.68;
                }
            }
            if(k==1){
                j++;
                i--;
            }
            else{
                i+=2;
            }
        }
		return energy;
	}
	
	private double[] calculateDNAEnergy(String mGene){
		double dH=0,dS=0,dG=0;
		int idx=0;
		boolean isInitial=true;
		boolean isTerminal=false;
		int numPairBases=0;
		int totalOverhang=0;
		int numOverhangs=0;
		int numInfBulgue=0;
		int numSupBulgue=0;
		for(int i=3; i<complex.size(); i+=3){
			//Overhang
			if(complex.get(i).equals(":")){
				idx+=Integer.parseInt(complex.get(i+1).toString());
				numOverhangs+=Integer.parseInt(complex.get(i+1).toString());
				double[] energy = fiveOverhangEnergy(mGene.charAt(idx-1), mGene.charAt(idx));
				dH+=energy[0]; dG+=energy[1]; dS += energy[1] * 1000 / 310.15;
				if(Integer.parseInt(complex.get(i+1).toString())>10){
					/*
					 * If Overhang is greater than 10, penalize it 
					 */
					dH+=Math.abs(-30.00*(Integer.parseInt(complex.get(i+1).toString())/10));
					dG+=Math.abs(-4.69*(Integer.parseInt(complex.get(i+1).toString())/10));
					dS+=Math.abs(-81.69*(Integer.parseInt(complex.get(i+1).toString())/10));
				}
				isInitial=false;
			}
			if(complex.get(i).equals("[")){
				if(i>=complex.size()-3){
					isTerminal=true;
				}
				double[] energy = cDuplexEnergy(mGene.substring(idx, idx + Integer.parseInt(complex.get(i+1).toString())), isInitial, isTerminal);
				idx+=Integer.parseInt(complex.get(i+1).toString());
				numPairBases+=Integer.parseInt(complex.get(i+1).toString());
				dH+=energy[0]; dS+=energy[1]; dG+=energy[2];
				isInitial=false;
			}
			if(complex.get(i).equals("(")){
				double energy=0;
				if(Integer.parseInt(complex.get(i+1).toString()) < 3){
					energy = internalSingleMiss(mGene.substring(idx - 1, idx + Integer.parseInt(complex.get(i+1).toString()) + 1));
				}
				else{
					energy=internalLoopGEnergy(Integer.parseInt(complex.get(i+1).toString()));
				}
				dG+=energy;
				dS+=energy*1000/310.15;
				idx+=Integer.parseInt(complex.get(i+1).toString());
				isInitial=false;
			}
			if(complex.get(i).equals("{") || complex.get(i).equals("<")){
				double energy=bulgesEnergy(Integer.parseInt(complex.get(i+1).toString()));
				dG+=energy;
				dS+=energy*1000/310.15;
				if(complex.get(i).equals("<")){
					numSupBulgue+=Integer.parseInt(complex.get(i+1).toString());
					idx+=Integer.parseInt(complex.get(i+1).toString());
				}
				else
					numInfBulgue += Integer.parseInt(complex.get(i+1).toString());
				isInitial=false;
			}
			if(complex.get(i).equals(";")){
				double[] energy = threeOverhangEnergy(mGene.charAt(idx-1), mGene.charAt(idx));
				idx+=Integer.parseInt(complex.get(i+1).toString());
				numOverhangs+=Integer.parseInt(complex.get(i+1).toString());
				dH+=energy[0]; dG+=energy[1]; dS += energy[1] * 1000 / 310.15;
				if(Integer.parseInt(complex.get(i+1).toString())>10){
					/*
					 * If Overhang is greater than 10, penalize it 
					 */
					dH+=Math.abs(-30.00*(Integer.parseInt(complex.get(i+1).toString())/10));
					dG+=Math.abs(-4.69*(Integer.parseInt(complex.get(i+1).toString())/10));
					dS+=Math.abs(-81.69*(Integer.parseInt(complex.get(i+1).toString())/10));
				}
			}
		}
		int shortestLen = 0;
		if(numInfBulgue > 0 && ( (numOverhangs + numSupBulgue) > numInfBulgue)){
			totalOverhang = numInfBulgue;
			shortestLen = Integer.parseInt(complex.get(1).toString())-Integer.parseInt(
					complex.get(0).toString()) + 1 - (numOverhangs + numSupBulgue) + totalOverhang;
		}
		else if(numInfBulgue > 0 && ( (numOverhangs + numSupBulgue) < numInfBulgue)){
			totalOverhang = numInfBulgue;
			shortestLen = Integer.parseInt(complex.get(1).toString())-Integer.parseInt(
					complex.get(0).toString()) + 1;
		}
		else{
			shortestLen = Integer.parseInt(complex.get(1).toString())-Integer.parseInt(
					complex.get(0).toString()) + 1 - (numOverhangs + numSupBulgue);
		}
		double Tm=dH*1000/(dS+R*Math.log(Ct/4))-273.15;
		//Calculate Kinetics
		this.numPairBases=numPairBases;
		this.totalOverHang=this.len - numPairBases;
		double[] kinetics = calculateKinetics(shortestLen, this.len, dG, Tm);
		double[] pEnergy= {dH,dS,dG,Tm,kinetics[0],kinetics[1]};//predicted Energy
		return pEnergy;
	}
	
	public double[] calculateHealthyRna_DxEnergy() {
        ArrayList<Object> hComplex = h_mRNAComplex();
        int idxSeq=0;
        int nBasePair = 0;
        boolean isInitial=true;
        String h_mRNA = this.h_mRNA.substring(
                Integer.parseInt(this.complex.get(0).toString()),
                Integer.parseInt(this.complex.get(1).toString()) + 1);
        double dH = 0, dS = 0, dG = 0;
        for(int i = 3; i < hComplex.size(); i+=3){
            if(hComplex.get(i).equals(":")){
                int nextNuc = Integer.parseInt(hComplex.get(i + 1).toString());
                nBasePair+=nextNuc;
                idxSeq += Integer.parseInt(hComplex.get(i + 1).toString());
                double[] energy = fiveOverhangRNAEnergy(h_mRNA.charAt(idxSeq - 1),
                        h_mRNA.charAt(idxSeq));
                dS += energy[1] * 1000 / 310.15;
                dG += energy[1];
                dH += energy[0];
            }
            else if(hComplex.get(i).equals("[")){
                int nextNuc = Integer.parseInt(hComplex.get(i+1).toString());
                double[] energy = cDuplexRNAEnergy(h_mRNA.substring(idxSeq, 
                        idxSeq + nextNuc), isInitial, true);
                idxSeq+=nextNuc;
                dH += energy[0];
                dS += energy[1];
                dG += energy[2];
                isInitial = false;
                nBasePair+=nextNuc;
                //if(nextNuc > nBasePair)
                    //nBasePair=nextNuc;
            }
            else if(hComplex.get(i).equals("(")){
                int nextNuc = Integer.parseInt(hComplex.get(i + 1).toString());
                double energy = 0;
                if (Integer.parseInt(hComplex.get(i + 1).toString()) < 3) {
                    energy = internalSingleMissRNA(h_mRNA.substring(idxSeq - 1, 
                            idxSeq + Integer.parseInt(hComplex.get(i + 1).toString()) + 1));
                } else {
                    energy = internalRNALoopGEnergy(Integer.parseInt(
                            hComplex.get(i + 1).toString()));
                }
                dG += energy;
                dS += energy * 1000 / 310.15;
                idxSeq += Integer.parseInt(hComplex.get(i + 1).toString());
                nBasePair+=nextNuc;
            }
            else if(hComplex.get(i).equals(";")){
                int nextNuc = Integer.parseInt(hComplex.get(i + 1).toString());
                nBasePair+=nextNuc;
                double[] energy = threeOverhangRNAEnergy(h_mRNA.charAt(idxSeq - 1),
                        h_mRNA.charAt(idxSeq));
                idxSeq += Integer.parseInt(hComplex.get(i + 1).toString());
                dS += energy[1] * 1000 / 310.15;
                dG += energy[1];
                dH += energy[0];
            }
        }
        double Tm = ((dH * 1000) / (dS + (R * Math.log(Ct / 4)))) - 273.15;
        if(dG>0.0)
            Tm=26;
        else if(Tm<26)
            Tm=26;
        double[] kinetics = calculateRNAKinetics(nBasePair,
                this.h_mRNA.length(), dG, Tm);
        //Predicted Energy
        double[] pEnergy = {dH, dS, dG, Tm, kinetics[0], kinetics[1]};
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
	private double[] calculateRNADNAEnergy(String mRNA, boolean isInitial){
		boolean isTerminal = true;
		double[] energy = cDuplexRNAEnergy(mRNA, isInitial, isTerminal);
		double Tm = energy[0] * 1000 / (energy[1] + R * Math.log(Ct / 4)) - 273.15;
		double[] kinetics = calculateKinetics(this.len, this.mGeneLen, energy[2], Tm);
		//double[] kinetics = calculateKinetics(this.len, this.totalOverHang, energy[2], Tm);
		double[] pEnergy= {energy[0],energy[1],energy[2],Tm,kinetics[0],kinetics[1]};//predicted Energy
		return pEnergy;
	}
	
	
	private double[] calculateKinetics(int shortestLen, int numPairBases, double dG, double Tm){
		//double incTemp=(Tm+273.15)-298.15; //Tm+273.15
		double incTemp = 37;
		double kF=(this.nRatio*Math.sqrt(shortestLen))/numPairBases;
		double kR=kF*Math.exp(dG/(this.R*incTemp));
		double[] kinectics = {kF, kR};
		return kinectics;
	}
	
	private double[] calculateRNAKinetics(int shortestLen, int numPairBases, double dG, double Tm) {
        //double incTemp = (Tm + 273.15) - 298.15; //Tm+273.15
        double incTemp = 37;
        //TODO: Revisar constante de nucleación en ARN
        double rnaNRatio = 1.3*Math.pow(10, 5);
        double kF = (this.nRatio * Math.sqrt(shortestLen)) / numPairBases;
        double kR = kF * Math.exp(dG / (this.R * incTemp));
        //double kOrderOne = Math.exp(dG/(this.R * incTemp));
        double[] kinectics = {kF, kR};
        return kinectics;
    }
    
    public String alignSequence(){
        StringBuilder alignedSequence = new StringBuilder();
        for(int i = Integer.parseInt(this.complex.get(0).toString()); 
                i <= Integer.parseInt(this.complex.get(1).toString()); i++){
            if(isRNA_DNAComplementary(this.mGene.charAt(i),
                    this.h_mRNA.charAt(i))){
                alignedSequence.append(this.h_mRNA.charAt(i));
            }
            else
                alignedSequence.append("-");
        }
        return alignedSequence.toString();
    }
    
    private ArrayList<Object> h_mRNAComplex(){
        ArrayList<Object> complexConfig = new ArrayList<Object>();
        int nContBase = 0;
        int nDisContBase = 0;
        boolean isInitial = true;
        boolean isTerminal = true;
        for(int i = Integer.parseInt(this.complex.get(0).toString()); 
                i <= Integer.parseInt(this.complex.get(1).toString()); i++){
            if(isRNA_DNAComplementary(this.mGene.charAt(i),
                    this.h_mRNA.charAt(i))){
                nContBase++;
                if (nDisContBase > 0){
                    if(isInitial){
                        complexConfig.add(":");
                        complexConfig.add(nDisContBase);
                        complexConfig.add(":");
                        isInitial = false;
                    }
                    else{
                        complexConfig.add("(");
                        complexConfig.add(nDisContBase);
                        complexConfig.add(")");
                    }
                }
                nDisContBase = 0;
                isTerminal = false;
            }
            else{
                if(nContBase > 0){
                    complexConfig.add("[");
                    complexConfig.add(nContBase);
                    complexConfig.add("]");
                    isInitial=false;
                }
                nDisContBase++;
                nContBase = 0;
                isTerminal = true;
            }
        }
        if(isTerminal){
            complexConfig.add(";");
            complexConfig.add(nDisContBase);
            complexConfig.add(";");
        }
        complexConfig.add(0, this.complex.get(0));
        complexConfig.add(1, this.complex.get(1));
        complexConfig.add(2, this.complex.get(2));
        //System.out.println(complexConfig.toString());
        return complexConfig;
    }
    
    private boolean isRNA_DNAComplementary(char a, char b){
        if(a == 'A' && b == 'U')
            return true;
        else if(a == 'U' && b == 'A')
            return true;
        else if(a == 'G' && b == 'C')
            return true;
        else if(a == 'C' && b == 'G')
            return true;
        else if(a == 'U' && b == 'A')
            return true;
        else if(a == 'A' && b == 'T')
            return true;
        else if(a == 'T' && b == 'A')
            return true;
        else
            return false;
    }
}
