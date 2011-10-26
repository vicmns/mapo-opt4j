package org.mapo;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import javax.swing.JFileChooser;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;


public class FileManager {
	
	  private static FileManager instance = null;
	   protected FileManager() {
	      // Exists only to defeat instantiation.
	   }
	   
	   public static FileManager getInstance() {
	      if(instance == null) {
	         instance = new FileManager();
	      }
	      return instance;
	   }
	
	public  String readStringFromFile(){
		String stringFromFile="";
		try{
			
			JFileChooser fileChooser = new JFileChooser();
			int returnOption = fileChooser.showOpenDialog(null);
			if(returnOption==JFileChooser.APPROVE_OPTION){
				File stringFile = fileChooser.getSelectedFile();
				FileInputStream fStream = new FileInputStream(stringFile);
				DataInputStream inStream = new DataInputStream(fStream);
				BufferedReader bReader = new BufferedReader(new InputStreamReader(inStream));
				//Read Line by line
				String strLine;
				while((strLine = bReader.readLine()) != null){
					strLine = strLine.replaceAll("\\s", "");
					strLine = strLine.replaceAll("-", "");
					stringFromFile += strLine;
				}
				return stringFromFile;
			}
		}
		catch(Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		return stringFromFile;
	}
}
