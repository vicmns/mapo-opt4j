package org.mapo;


import com.google.inject.Inject;

public class ComplexProblem {

	protected int idxS;
	  protected int idxE;
	  protected String Gene;
	  protected String AB;
	  protected int mGeneLen;
	  protected String mRNA;
	  protected String h_mRNA;

	  @Inject
	  public ComplexProblem()
	  {
		//TODO: Delete path from readFromFile to let it choose
	    String gene = readFromFile("/home/vicmns/Sequence/CFTR-cDNA-1399-2004");
		//String gene = readFromFile(System.getProperty("user.home") + "\\Documents\\CFTR-cDNA-1399-2004");
		//String gene = readFromFile();
	    if (gene == "")
	      return;
	    String ab = readFromFile("/home/vicmns/Sequence/CFTR-cDNA-1399-2004-dF508(c.1521_1523delCTT)");
	    //String ab = readFromFile();
	    //String ab = readFromFile(System.getProperty("user.home") + "\\Documents\\CFTR-cDNA-1399-2004-dF508(c.1521_1523delCTT)");
	    if (ab == "")
	      return;
	    this.Gene = gene;
	    this.AB = ab;
	    this.mGeneLen = ab.length();
	    this.mRNA = cDNAToRNA(this.AB);
	    this.h_mRNA = cDNAToRNA(this.Gene);
	    int[] idx = findMutation(this.Gene, this.AB);
	    this.idxS = idx[0];
	    this.idxE = idx[1];
	  }
	  
	  public String getHmRNA() {
		  return this.h_mRNA;
	  }

	  public int getMGeneLen() {
	    return this.mGeneLen;
	  }

	  public int getIdxS() {
	    return this.idxS;
	  }

	  public int getIdxE() {
	    return this.idxE;
	  }

	  public int getGeneLength() {
	    return this.Gene.length();
	  }

	  public int getABLength() {
	    return this.AB.length();
	  }

	  public String getGene() {
	    return this.Gene;
	  }

	  public String getAB() {
	    return this.AB;
	  }

	  public String getMRNA() {
	    return this.mRNA;
	  }

	  private String readFromFile() {
	    FileManager fileMan = new FileManager();
	    return fileMan.readStringFromFile();
	  }
	  private String readFromFile(String srcFile) {
	    FileManager fileMan = new FileManager();
	    return fileMan.readStringFromFile(srcFile);
	  }

	  private int[] findMutation(String a, String b) {
	    int firstPos = 0;
	    int[] idx = { -1, -1 };
	    int len = 0;
	    boolean hasMutation = true;
	    int i = 0;
	    if (a.length() > b.length())
	      len = b.length();
	    else
	      len = a.length();
	    while ((hasMutation) && (i < len)) {
	      if (a.charAt(i) != b.charAt(i)) {
	        if (idx[0] < 0)
	          idx[0] = i;
	        if (a.charAt(i + 1) == b.charAt(i + 1)) {
	          hasMutation = false;
	        }
	        firstPos++;
	      }
	      i++;
	    }
	    idx[1] = (idx[0] + firstPos - 1);
	    return idx;
	  }

	  private String cDNAToRNA(String sequence) {
	    String RNAString = "";
	    for (int i = 0; i < sequence.length(); i++) {
	      switch (sequence.charAt(i)) {
	      case 'A':
	        RNAString = RNAString + "U";
	        break;
	      case 'T':
	        RNAString = RNAString + "A";
	        break;
	      case 'C':
	        RNAString = RNAString + "G";
	        break;
	      case 'G':
	        RNAString = RNAString + "C";
	      }
	    }

	    return RNAString;
	  }
}
