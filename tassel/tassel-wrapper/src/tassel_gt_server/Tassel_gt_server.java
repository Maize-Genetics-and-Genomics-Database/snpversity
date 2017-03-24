/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tassel_gt_server;

/**
 *
 * @author qs24
 * 
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.log4j.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.JSONValue;
import java.io.FileWriter;
import java.io.IOException;

/* Boon Stuff */
import org.boon.Boon;
import org.boon.IO;
import org.boon.json.JsonFactory;
import org.boon.json.JsonParserAndMapper;
import org.boon.json.JsonParserFactory;
import org.boon.json.ObjectMapper;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

/*Boon End */
import com.google.gson.stream.JsonWriter;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;




public class Tassel_gt_server {

	private static final Logger myLogger = Logger.getLogger(Tassel_gt_server.class);
	
    /**
     * @param args the command line arguments
     * create_result_file T or F
     * 
     */
    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.print("No method is specified!\n");
            printUsage("main");
        }
        String method = args[0];
        String[] newargs = new String[args.length-1];
        System.arraycopy(args, 1, newargs, 0, args.length-1);
        if (method.equals("slice")) 
        {
            String results = slice(newargs);
            System.out.print(results);
        }
        else if (method.equals("sliceNuke")) 
        {
            String results = sliceNuke(newargs);
            System.out.print(results);
        }
        else if (method.equals("sliceRangeNuke")) 
        {
            String results = sliceRangeNuke(newargs);
            System.out.print(results);
        }
        else if (method.equals("dbinfo")) 
        {
            String results = dbinfo(newargs);
            System.out.print(results);
        }
        else if (method.equals("get_taxa_list"))
        {
            String[] results = get_taxa_list(newargs);
            for (String t:results)
            {
                System.out.println(t);
            }
        }
        else if (method.equals("chr_marker_info"))
        {
            String[] results = chr_marker_info(newargs);
            for (String t:results)
            {
                System.out.println(t);
            }
        }
        else if (method.equals("sliceSplitFile"))
        {
            String results = sliceSplitFile(newargs);
            System.out.print(results);
        }
        else
        {
           System.out.print("The method '"+ method + "' is not recognized!\n");
           printUsage("main");
        }
    }
    
    private static String sliceNuke(String[] args)
    {
    	
    	//get parameters
        ArgsEngine  myArgsEngine = new ArgsEngine();
        myArgsEngine.add("-sf", "--sourcefile", true);
        myArgsEngine.add("-st", "--sourcefile-type", true);
        myArgsEngine.add("-df", "--destinationfile", true);
        myArgsEngine.add("-dt", "--destinationfile-type", true);
        myArgsEngine.add("-tf", "--taxa-file", true);
        myArgsEngine.add("-tl", "--taxa-list", true);
        myArgsEngine.add("-ch", "--chromosome", true);
        myArgsEngine.add("-start", "--chr-start", true);
        myArgsEngine.add("-end", "--chr-end", true);
        myArgsEngine.add("-bf", "--build-file", false);
        myArgsEngine.parse(args);
     
       boolean buildfile = false;
       String source_file = null;
       String source_file_type = null;
       String dest_file = null;
       String dest_file_type=null;
       String TaxaListFile=null;
       ArrayList<String> TaxaArrayList = new ArrayList<String>();
       String ChromosomeStr=null;
       int StartPhysicalPosition = 0;
       int EndPhysicalPosition = 0;

       if (myArgsEngine.getBoolean("-bf")) 
       {
            buildfile = true;      
       }
       
      if (myArgsEngine.getBoolean("-sf")) {
          source_file = myArgsEngine.getString("-sf");
      } else {
          printUsage("slice");
          throw new IllegalArgumentException("Please specify a source file (option -sf).");
      }
      
      if (myArgsEngine.getBoolean("-st")) {
          source_file_type = myArgsEngine.getString("-st");
      } else {
          printUsage("slice");
          throw new IllegalArgumentException("Please specify a source file type (option -st).");
      }
      
      if (buildfile)
      {
          if (myArgsEngine.getBoolean("-df")) {
              dest_file = myArgsEngine.getString("-df");
          } else {
              printUsage("slice");
              throw new IllegalArgumentException("Please specify a destination file (option -df).");
          }

          if (myArgsEngine.getBoolean("-dt")) {
              dest_file_type = myArgsEngine.getString("-dt");
          } else {
              printUsage("slice");
              throw new IllegalArgumentException("Please specify a destination file type (option -dt).");
          }
      }
      
      if (myArgsEngine.getBoolean("-tf")) 
      {
          TaxaListFile = myArgsEngine.getString("-tf");
          File outDirectory = new File(TaxaListFile);
          if (!outDirectory.isFile()) {
              printUsage("slice");
              throw new IllegalArgumentException("The taxa file you supplied (option -tf) is not a file: " + TaxaListFile);
          }
          //verify and create sub-taxalist
          //create taxa filtered genotype table
         
         try {
              BufferedReader br = new BufferedReader(new FileReader(TaxaListFile), 65536);

              String temp;
              int currLine = 0;
              while (((temp = br.readLine()) != null)) {
                  if (!temp.trim().isEmpty())
                  {
                      TaxaArrayList.add(temp.trim());
                      currLine++;
                  }

              }
          } catch (Exception e) {
              System.out.println("Couldn't open taxa file to read taxa list: " + e);
          }
      } 
      else if (myArgsEngine.getBoolean("-tl"))
      {
          String TaxaListStr = myArgsEngine.getString("-tf");
          if (TaxaListStr.equalsIgnoreCase("all"))
          {
              TaxaArrayList.add("all");
          }
          else
          {
              String[] wordList = TaxaListStr.split(";");
              TaxaArrayList.addAll(Arrays.asList(wordList));  
          }
      }
      else {
          printUsage("slice");
          throw new IllegalArgumentException("Please specify the file with taxa list (option -tf).");
      }
                             
      if (myArgsEngine.getBoolean("-ch")) {
          ChromosomeStr = myArgsEngine.getString("-ch");
      } else {
          if ((myArgsEngine.getBoolean("-start")) || (myArgsEngine.getBoolean("-end")))
          {
              printUsage("slice");
              throw new IllegalArgumentException("Please specify a chromosome name (option -ch).");
          }
          else
          {
              ChromosomeStr = "all";
          }
      }
      
      if (myArgsEngine.getBoolean("-start")) {
          StartPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-start"));
      } 
      
      if (myArgsEngine.getBoolean("-end")) {
          EndPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-end"));
      }

       //create a genotypetable from source
      GenotypeTable source_gt_table = null;
      if (source_file_type.equalsIgnoreCase("hdf5"))
      {
          source_gt_table = GenotypeTableBuilder.getInstance(source_file);
      }
      else if (source_file_type.equalsIgnoreCase("hapmap"))
      {
          source_gt_table = ImportUtils.readFromHapmap(source_file);
      }
      else if (source_file_type.equalsIgnoreCase("vcf"))
      {
          source_gt_table = ImportUtils.readFromVCF(source_file, null, false);
      }

      

      
      GenotypeTable gt_taxa_filtered = null ; 
      
              
      
      //if TaxaListString=all or not specified, all taxa included
      if ((TaxaArrayList.get(0).equalsIgnoreCase("all")) || (TaxaArrayList.size()==0))
      {
          gt_taxa_filtered = source_gt_table;
      }
      //if TaxaListString specified, create subset
      else
      {
          TaxaList source_taxa_list = source_gt_table.taxa();
          TaxaList dest_taxa_list=null;
          TaxaListBuilder dest_taxa_list_builder = new TaxaListBuilder();

          ArrayList<String> unknow_taxa = new ArrayList<String>();
          for (int i =0; i<TaxaArrayList.size(); i++)
          {
             if (TaxaArrayList.get(i).matches(".*\\w.*"))
             {
                  int itaxa_index = source_taxa_list.indexOf(TaxaArrayList.get(i));
                  if (itaxa_index>=0)
                  {
                      dest_taxa_list_builder.add(source_taxa_list.get(itaxa_index));
                  }
                  else
                  {
                      unknow_taxa.add(TaxaArrayList.get(i));
                  }
             }
          }
          dest_taxa_list = dest_taxa_list_builder.build();
          if (unknow_taxa.size()>0)
          {
              String errmsg = "-1 -1\nUnknown taxa names:\n";
              for (int i=0; i<unknow_taxa.size(); i++)
              {
                  errmsg += unknow_taxa.get(i) + "\n";
              }
              return errmsg;
          }
          gt_taxa_filtered = FilterGenotypeTable.getInstance(source_gt_table, dest_taxa_list);
      }
      //create sites filtered genotype table
      GenotypeTable gt_taxa_sites_filtered = null;
      
      //create position filtered genotype table
      if (ChromosomeStr.toLowerCase().equalsIgnoreCase("all"))
      {
          gt_taxa_sites_filtered = gt_taxa_filtered;
      }
      else
      {
          Chromosome ch = source_gt_table.chromosome(ChromosomeStr);
          if (ch == null)
          {
              return "-1 -1\n" + "The chromosome name '" + ChromosomeStr +  "' is not recognized\n";
          }
          if ((StartPhysicalPosition==0) || (EndPhysicalPosition==0))
          {
              gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch);
          }
          else
          {
              gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch, StartPhysicalPosition, EndPhysicalPosition);
          }
          
      }
      	BufferedWriter writer = null;
      	try {
      		
	      		File outFile = new File(dest_file + ".txt");
	      		writer = new BufferedWriter(new FileWriter(outFile));
	            System.out.println(outFile.getCanonicalPath());
  	
	      		writer.write("Selected Taxa : " + gt_taxa_sites_filtered.taxaName(0));
	      		writer.write("\nSelected Chromosome : " + gt_taxa_sites_filtered.chromosome(ChromosomeStr).toString());
	      		writer.write("\nSelected Position : " + StartPhysicalPosition);
	      	      	
	      	// Get the Unique site number for the specific physical position and specific Chromosome.
	     	int site = source_gt_table.siteOfPhysicalPosition(StartPhysicalPosition, gt_taxa_sites_filtered.chromosome(ChromosomeStr));
	     	
	     	writer.write("\nSite Value for the selected Chromosome & Position : " + site);
	     	writer.write("\nGenotype Present at the selected Chromosome & Position in the selected Taxa : " + gt_taxa_sites_filtered.genotypeAsString(0,0));
	     	writer.write("\nGenotype byte value : " + gt_taxa_sites_filtered.genotype(0,0));
	      	
	     	byte[] allGenotypeDiploidValues = source_gt_table.genotypeAllTaxa(site);
	     	byte genotypeDiploidValueAtSpecificPosition = gt_taxa_sites_filtered.genotype(0,0);
	      	BitSet matchedTaxa = GenotypeTableUtils.calcBitPresenceOfDiploidValueFromGenotype(allGenotypeDiploidValues, genotypeDiploidValueAtSpecificPosition);
	      	int[] matchedPositions = matchedTaxa.getIndicesOfSetBits();
	      	writer.write("\nNumber of matches : " + matchedPositions.length);
	      	writer.write("\nList of all Taxa containing the same Genotype at the same Chromosome & Position : ");

	      	for (int i : matchedPositions)
		      	{
	       		writer.write("\n" + source_gt_table.taxaName(i));
		      	}
      	} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
      	finally
      	{
      		try {
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
      	}
      	System.out.println("Done");
      return "";
    	
    }
    
    private static String sliceRangeNuke(String[] args)
    {
    	//get parameters
        ArgsEngine  myArgsEngine = new ArgsEngine();
        myArgsEngine.add("-sf", "--sourcefile", true);
        myArgsEngine.add("-st", "--sourcefile-type", true);
        myArgsEngine.add("-df", "--destinationfile", true);
        myArgsEngine.add("-dt", "--destinationfile-type", true);
        myArgsEngine.add("-tf", "--taxa-file", true);
        myArgsEngine.add("-tl", "--taxa-list", true);
        myArgsEngine.add("-ch", "--chromosome", true);
        myArgsEngine.add("-start", "--chr-start", true);
        myArgsEngine.add("-end", "--chr-end", true);
        myArgsEngine.add("-bf", "--build-file", false);
        myArgsEngine.parse(args);
     
       boolean buildfile = false;
       String source_file = null;
       String source_file_type = null;
       String dest_file = null;
       String dest_file_type=null;
       String TaxaListFile=null;
       ArrayList<String> TaxaArrayList = new ArrayList<String>();
       String ChromosomeStr=null;
       int StartPhysicalPosition = 0;
       int EndPhysicalPosition = 0;

       if (myArgsEngine.getBoolean("-bf")) 
       {
            buildfile = true;      
       }
       
      if (myArgsEngine.getBoolean("-sf")) {
          source_file = myArgsEngine.getString("-sf");
      } else {
          printUsage("slice");
          throw new IllegalArgumentException("Please specify a source file (option -sf).");
      }
      
      if (myArgsEngine.getBoolean("-st")) {
          source_file_type = myArgsEngine.getString("-st");
      } else {
          printUsage("slice");
          throw new IllegalArgumentException("Please specify a source file type (option -st).");
      }
      
      if (buildfile)
      {
          if (myArgsEngine.getBoolean("-df")) {
              dest_file = myArgsEngine.getString("-df");
          } else {
              printUsage("slice");
              throw new IllegalArgumentException("Please specify a destination file (option -df).");
          }

          if (myArgsEngine.getBoolean("-dt")) {
              dest_file_type = myArgsEngine.getString("-dt");
          } else {
              printUsage("slice");
              throw new IllegalArgumentException("Please specify a destination file type (option -dt).");
          }
      }
      
      if (myArgsEngine.getBoolean("-tf")) 
      {
          TaxaListFile = myArgsEngine.getString("-tf");
          File outDirectory = new File(TaxaListFile);
          if (!outDirectory.isFile()) {
              printUsage("slice");
              throw new IllegalArgumentException("The taxa file you supplied (option -tf) is not a file: " + TaxaListFile);
          }
          //verify and create sub-taxalist
          //create taxa filtered genotype table
         
         try {
              BufferedReader br = new BufferedReader(new FileReader(TaxaListFile), 65536);

              String temp;
              int currLine = 0;
              while (((temp = br.readLine()) != null)) {
                  if (!temp.trim().isEmpty())
                  {
                      TaxaArrayList.add(temp.trim());
                      currLine++;
                  }

              }
          } catch (Exception e) {
              System.out.println("Couldn't open taxa file to read taxa list: " + e);
          }
      } 
      else if (myArgsEngine.getBoolean("-tl"))
      {
          String TaxaListStr = myArgsEngine.getString("-tf");
          if (TaxaListStr.equalsIgnoreCase("all"))
          {
              TaxaArrayList.add("all");
          }
          else
          {
              String[] wordList = TaxaListStr.split(";");
              TaxaArrayList.addAll(Arrays.asList(wordList));  
          }
      }
      else {
          printUsage("slice");
          throw new IllegalArgumentException("Please specify the file with taxa list (option -tf).");
      }
                             
      if (myArgsEngine.getBoolean("-ch")) {
          ChromosomeStr = myArgsEngine.getString("-ch");
      } else {
          if ((myArgsEngine.getBoolean("-start")) || (myArgsEngine.getBoolean("-end")))
          {
              printUsage("slice");
              throw new IllegalArgumentException("Please specify a chromosome name (option -ch).");
          }
          else
          {
              ChromosomeStr = "all";
          }
      }
      
      if (myArgsEngine.getBoolean("-start")) {
          StartPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-start"));
      } 
      
      if (myArgsEngine.getBoolean("-end")) {
          EndPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-end"));
      }

       //create a genotypetable from source
      GenotypeTable source_gt_table = null;
      if (source_file_type.equalsIgnoreCase("hdf5"))
      {
          source_gt_table = GenotypeTableBuilder.getInstance(source_file);
      }
      else if (source_file_type.equalsIgnoreCase("hapmap"))
      {
          source_gt_table = ImportUtils.readFromHapmap(source_file);
      }
      else if (source_file_type.equalsIgnoreCase("vcf"))
      {
          source_gt_table = ImportUtils.readFromVCF(source_file, null, false);
      }

      

      
      GenotypeTable gt_taxa_filtered = null ; 
      
              
      
      //if TaxaListString=all or not specified, all taxa included
      if ((TaxaArrayList.get(0).equalsIgnoreCase("all")) || (TaxaArrayList.size()==0))
      {
          gt_taxa_filtered = source_gt_table;
      }
      //if TaxaListString specified, create subset
      else
      {
          TaxaList source_taxa_list = source_gt_table.taxa();
          TaxaList dest_taxa_list=null;
          TaxaListBuilder dest_taxa_list_builder = new TaxaListBuilder();

          ArrayList<String> unknow_taxa = new ArrayList<String>();
          for (int i =0; i<TaxaArrayList.size(); i++)
          {
             if (TaxaArrayList.get(i).matches(".*\\w.*"))
             {
                  int itaxa_index = source_taxa_list.indexOf(TaxaArrayList.get(i));
                  if (itaxa_index>=0)
                  {
                      dest_taxa_list_builder.add(source_taxa_list.get(itaxa_index));
                  }
                  else
                  {
                      unknow_taxa.add(TaxaArrayList.get(i));
                  }
             }
          }
          dest_taxa_list = dest_taxa_list_builder.build();
          if (unknow_taxa.size()>0)
          {
              String errmsg = "-1 -1\nUnknown taxa names:\n";
              for (int i=0; i<unknow_taxa.size(); i++)
              {
                  errmsg += unknow_taxa.get(i) + "\n";
              }
              return errmsg;
          }
          gt_taxa_filtered = FilterGenotypeTable.getInstance(source_gt_table, dest_taxa_list);
      }
      //create sites filtered genotype table
      GenotypeTable gt_taxa_sites_filtered = null;
      
      //create position filtered genotype table
      if (ChromosomeStr.toLowerCase().equalsIgnoreCase("all"))
      {
          gt_taxa_sites_filtered = gt_taxa_filtered;
      }
      else
      {
          Chromosome ch = source_gt_table.chromosome(ChromosomeStr);
          if (ch == null)
          {
              return "-1 -1\n" + "The chromosome name '" + ChromosomeStr +  "' is not recognized\n";
          }
          if ((StartPhysicalPosition==0) || (EndPhysicalPosition==0))
          {
              gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch);
          }
          else
          {
              gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch, StartPhysicalPosition, EndPhysicalPosition);
          }
          
      }
      	BufferedWriter writer = null;
      	try {
      		
	      		File outFile = new File(dest_file + ".txt");
	      		writer = new BufferedWriter(new FileWriter(outFile));
	            System.out.println(outFile.getCanonicalPath());
  	
	      		writer.write("Selected Taxa : " + gt_taxa_sites_filtered.taxaName(0));
	      		writer.write("\nSelected Chromosome : " + gt_taxa_sites_filtered.chromosome(ChromosomeStr).toString());
	      		writer.write("\nSelected Positions : ");
	      		ArrayList sites = new ArrayList<Integer>();
	      		int c=0;
	      		for (int i : gt_taxa_sites_filtered.physicalPositions())
	      		{
	      			sites.add(source_gt_table.siteOfPhysicalPosition(i, gt_taxa_sites_filtered.chromosome(ChromosomeStr)));
	      			writer.write("\n" + i + " - Site Value - " + sites.get(c) +" - Genotype present - " + gt_taxa_sites_filtered.genotypeAsString(0,c) + " - Byte Value - " + gt_taxa_sites_filtered.genotype(0,c));
	      			c++;
	      		}
	      		    
	      		ArrayList<BitSet> allMatchs = new ArrayList<BitSet>();

	      		for (int i = 0; i<c; i++)
	      		{
	      			byte genotypeDiploidValueAtSpecificPosition = gt_taxa_sites_filtered.genotype(0,i);
	      			if (genotypeDiploidValueAtSpecificPosition == -1) break;
	      			byte[] allGenotypeDiploidValues = source_gt_table.genotypeAllTaxa((int) sites.get(i));
	      			BitSet matchedTaxa = GenotypeTableUtils.calcBitPresenceOfDiploidValueFromGenotype(allGenotypeDiploidValues, genotypeDiploidValueAtSpecificPosition);
	      			allMatchs.add(matchedTaxa);
	      		}
	      		
	      		BitSet globalMatch = allMatchs.get(0);
	      		for (int i=1; i<allMatchs.size(); i++)
	      		{
	      			globalMatch.and(allMatchs.get(i));
	      		}
	      		
	      	int[] matchedPositions = globalMatch.getIndicesOfSetBits();
	      	writer.write("\nNumber of matches : " + matchedPositions.length);
	      	writer.write("\nList of all Taxa containing the same Genotype at the same Chromosome & Position : ");

	      	for (int i : matchedPositions)
		      	{
	       		writer.write("\n" + source_gt_table.taxaName(i));
		      	}
      	} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
      	finally
      	{
      		try {
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
      	}
      	System.out.println("Done");
      return "";
    	
    }
    
    private static String dbinfo(String[] args)
    {
          //get parameters
          ArgsEngine  myArgsEngine = new ArgsEngine();
          myArgsEngine.add("-sf", "--sourcefile", true);
          myArgsEngine.add("-st", "--sourcefile-type", true);
          myArgsEngine.add("-dbname", "--dbname", true);
          myArgsEngine.add("-o", "--out", true);
          myArgsEngine.parse(args);
          String source_file ="";
          String source_file_type = "";
          String dbname = "";
          String outdir = "";
         if (myArgsEngine.getBoolean("-sf")) {
            source_file = myArgsEngine.getString("-sf");
          } else {
            printUsage("dbinfo");
            throw new IllegalArgumentException("Please specify a source file (option -sf).");
        }
        
        if (myArgsEngine.getBoolean("-st")) {
            source_file_type = myArgsEngine.getString("-st");
        } else {
            printUsage("dbinfo");
            throw new IllegalArgumentException("Please specify a source file type (option -st).");
        }
        
        if (myArgsEngine.getBoolean("-dbname")) {
            dbname = myArgsEngine.getString("-dbname");
        } else {
            printUsage("dbinfo");
            throw new IllegalArgumentException("Please specify a database name (option -dbname ).");
        }
        
        if (myArgsEngine.getBoolean("-o")) {
            outdir = myArgsEngine.getString("-o");
        } else {
            outdir = ".";
        }
                
        //create a genotypetable from source
        GenotypeTable source_gt_table = null;
        if (source_file_type.equalsIgnoreCase("hdf5"))
        {
            source_gt_table = GenotypeTableBuilder.getInstance(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("hapmap"))
        {
            source_gt_table = ImportUtils.readFromHapmap(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("vcf"))
        {
            source_gt_table = ImportUtils.readFromVCF(source_file, null, false);
        }
        
        
        BufferedWriter bw = null;
        //write the taxa file
        try {        
            bw = Utils.getBufferedWriter(outdir + "/" + dbname + ".taxainfo");
            TaxaList taxalist  = source_gt_table.taxa();
            for (int i=0; i<taxalist.size(); i++)
            {
                Taxon t = taxalist.get(i);
                Map.Entry<String,String>[] tannotations = t.getAnnotation().getAllAnnotationEntries(); // t.getAllAnnotationEntries();
                if (tannotations.length>0)
                {
                    for (Map.Entry<String,String> ta:tannotations)
                    {
                          bw.write(dbname + "\t" + t.getName() + "\t" + ta.getKey() + "\t" + ta.getValue() + "\n");
                    }
                }
                else
                {
                    bw.write(dbname + "\t" + t.getName() + "\t" + "" + "\t" + "" + "\n");
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing dbinfo taxainfo file: " + outdir + "/" + dbname + ".taxainfo " + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        //write chromosome file
        try {        
            bw = Utils.getBufferedWriter(outdir + "/" + dbname + ".chrinfo");
            Chromosome[] chrlist  = source_gt_table.chromosomes();
            for (int i=0; i<chrlist.length; i++)
            {
                Chromosome c = chrlist[i];
                int sitecount = source_gt_table.chromosomeSiteCount(c);
                //int[] firstlast = source_gt_table.firstLastSiteOfChromosome(c);
                bw.write("\t" + dbname + "\t" + c.getName() + "\t" + c.getLength() + "\t" + sitecount  + "\n");
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing chrinfo file: " + outdir + "/" + dbname + ".chrinfo " + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        //write dbinfo
        try {        
            bw = Utils.getBufferedWriter(outdir + "/" + dbname + ".dbinfo");
            String genomeversion = source_gt_table.genomeVersion();
            int hasdepth =0;
            if (source_gt_table.hasDepth())
            {
                hasdepth=1;
            }
            String date_created = "";
            String tassel_version = "";
            String description = "DB description unavailable.";
            bw.write("\t" + dbname + "\t" + description + "\t" + genomeversion + "\t" + source_file + "\t" + source_file_type + "\t" + date_created + "\t" + hasdepth + "\t" + tassel_version + "\n" );
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing dbinfo dbinfo file: " + outdir + "/" + dbname + "dbinfo " + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return "done";
        
    }
    private static String sliceSplitFile(String[] args)
    {
          //get parameters
          ArgsEngine  myArgsEngine = new ArgsEngine();
          myArgsEngine.add("-sf", "--sourcefile", true);
          myArgsEngine.add("-st", "--sourcefile-type", true);
          myArgsEngine.add("-df", "--destinationfile", true);
          myArgsEngine.add("-dt", "--destinationfile-type", true);
          myArgsEngine.add("-tf", "--taxa-file", true);
          myArgsEngine.add("-tl", "--taxa-list", true);
          myArgsEngine.add("-ch", "--chromosome", true);
          myArgsEngine.add("-start", "--chr-start", true);
          myArgsEngine.add("-end", "--chr-end", true);
          myArgsEngine.add("-bf", "--build-file", false);
          myArgsEngine.add("-split", "--split-by-sites", true);
          myArgsEngine.parse(args);
       
         boolean buildfile = false;
         String source_file = null;
         String source_file_type = null;
         String dest_file = null;
         String dest_file_type=null;
         String TaxaListFile=null;
         ArrayList<String> TaxaArrayList = new ArrayList<String>();
         String ChromosomeStr=null;
         int StartPhysicalPosition = 0;
         int EndPhysicalPosition = 0;
         int splitAfterSite = 0;

         if (myArgsEngine.getBoolean("-bf")) 
         {
              buildfile = true;      
         }
         
        if (myArgsEngine.getBoolean("-sf")) {
            source_file = myArgsEngine.getString("-sf");
        } else {
            printUsage("slice");
            throw new IllegalArgumentException("Please specify a source file (option -sf).");
        }
        
        if (myArgsEngine.getBoolean("-st")) {
            source_file_type = myArgsEngine.getString("-st");
        } else {
            printUsage("slice");
            throw new IllegalArgumentException("Please specify a source file type (option -st).");
        }
        
        if (buildfile)
        {
            if (myArgsEngine.getBoolean("-df")) {
                dest_file = myArgsEngine.getString("-df");
            } else {
                printUsage("slice");
                throw new IllegalArgumentException("Please specify a destination file (option -df).");
            }

            if (myArgsEngine.getBoolean("-dt")) {
                dest_file_type = myArgsEngine.getString("-dt");
            } else {
                printUsage("slice");
                throw new IllegalArgumentException("Please specify a destination file type (option -dt).");
            }
        }
        
        if (myArgsEngine.getBoolean("-tf")) 
        {
            TaxaListFile = myArgsEngine.getString("-tf");
            File outDirectory = new File(TaxaListFile);
            if (!outDirectory.isFile()) {
                printUsage("slice");
                throw new IllegalArgumentException("The taxa file you supplied (option -tf) is not a file: " + TaxaListFile);
            }
            //verify and create sub-taxalist
            //create taxa filtered genotype table
           
           try {
                BufferedReader br = new BufferedReader(new FileReader(TaxaListFile), 65536);

                String temp;
                int currLine = 0;
                while (((temp = br.readLine()) != null)) {
                    if (!temp.trim().isEmpty())
                    {
                        TaxaArrayList.add(temp.trim());
                        currLine++;
                    }

                }
            } catch (Exception e) {
                System.out.println("Couldn't open taxa file to read taxa list: " + e);
            }
        } 
        else if (myArgsEngine.getBoolean("-tl"))
        {
            String TaxaListStr = myArgsEngine.getString("-tf");
            if (TaxaListStr.equalsIgnoreCase("all"))
            {
                TaxaArrayList.add("all");
            }
            else
            {
                String[] wordList = TaxaListStr.split(";");
                TaxaArrayList.addAll(Arrays.asList(wordList));  
            }
        }
        else {
            printUsage("slice");
            throw new IllegalArgumentException("Please specify the file with taxa list (option -tf).");
        }
                               
        if (myArgsEngine.getBoolean("-ch")) {
            ChromosomeStr = myArgsEngine.getString("-ch");
        } else {
            if ((myArgsEngine.getBoolean("-start")) || (myArgsEngine.getBoolean("-end")))
            {
                printUsage("slice");
                throw new IllegalArgumentException("Please specify a chromosome name (option -ch).");
            }
            else
            {
                ChromosomeStr = "all";
            }
        }
        
        if (myArgsEngine.getBoolean("-start")) {
            StartPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-start"));
        } 
        
        if (myArgsEngine.getBoolean("-end")) {
            EndPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-end"));
        }
        if (myArgsEngine.getBoolean("-split")) {
            splitAfterSite = Integer.parseInt(myArgsEngine.getString("-split"));
        }

         //create a genotypetable from source
        GenotypeTable source_gt_table = null;
        if (source_file_type.equalsIgnoreCase("hdf5"))
        {
            source_gt_table = GenotypeTableBuilder.getInstance(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("hapmap"))
        {
            source_gt_table = ImportUtils.readFromHapmap(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("vcf"))
        {
            source_gt_table = ImportUtils.readFromVCF(source_file, null, false);
        }

        

        
        GenotypeTable gt_taxa_filtered = null ; 
        
                
        
        //if TaxaListString=all or not specified, all taxa included
        if ((TaxaArrayList.get(0).equalsIgnoreCase("all")) || (TaxaArrayList.size()==0))
        {
            gt_taxa_filtered = source_gt_table;
        }
        //if TaxaListString specified, create subset
        else
        {
            TaxaList source_taxa_list = source_gt_table.taxa();
            TaxaList dest_taxa_list=null;
            TaxaListBuilder dest_taxa_list_builder = new TaxaListBuilder();

            ArrayList<String> unknow_taxa = new ArrayList<String>();
            for (int i =0; i<TaxaArrayList.size(); i++)
            {
               if (TaxaArrayList.get(i).matches(".*\\w.*"))
               {
                    int itaxa_index = source_taxa_list.indexOf(TaxaArrayList.get(i));
                    if (itaxa_index>=0)
                    {
                        dest_taxa_list_builder.add(source_taxa_list.get(itaxa_index));
                    }
                    else
                    {
                        unknow_taxa.add(TaxaArrayList.get(i));
                    }
               }
            }
            dest_taxa_list = dest_taxa_list_builder.build();
            if (unknow_taxa.size()>0)
            {
                String errmsg = "-1 -1\nUnknown taxa names:\n";
                for (int i=0; i<unknow_taxa.size(); i++)
                {
                    errmsg += unknow_taxa.get(i) + "\n";
                }
                return errmsg;
            }
            gt_taxa_filtered = FilterGenotypeTable.getInstance(source_gt_table, dest_taxa_list);
        }
        
        
        //create sites filtered genotype table
        GenotypeTable gt_taxa_sites_filtered = null;
        
        //create position filtered genotype table
        if (ChromosomeStr.toLowerCase().equalsIgnoreCase("all"))
        {
            gt_taxa_sites_filtered = gt_taxa_filtered;
        }
        else
        {
            Chromosome ch = source_gt_table.chromosome(ChromosomeStr);
            if (ch == null)
            {
                return "-1 -1\n" + "The chromosome name '" + ChromosomeStr +  "' is not recognized\n";
            }
            if ((StartPhysicalPosition==0) || (EndPhysicalPosition==0))
            {
                gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch);
            }
            else
            {
                gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch, StartPhysicalPosition, EndPhysicalPosition);
            }
            
        }
        String return_Values = "";
        if (gt_taxa_sites_filtered == null) 
        {
            return "-1 -1\nNo data after filtering!\n";
        }
        else
        {
            return_Values= gt_taxa_sites_filtered.taxa().size() + " " + gt_taxa_sites_filtered.positions().size();
        }
        if (buildfile)
        {
            if (dest_file_type.equalsIgnoreCase("hdf5"))
            {
                ExportUtils.writeGenotypeHDF5(gt_taxa_sites_filtered, dest_file, source_gt_table.hasDepth());
            }
            else if (dest_file_type.equalsIgnoreCase("hapmap"))
            {
                ExportUtils.writeToHapmap(gt_taxa_sites_filtered, dest_file);
            }
            else if (dest_file_type.equalsIgnoreCase("vcf"))
            {
                ExportUtils.writeToVCF(gt_taxa_sites_filtered, dest_file, source_gt_table.hasDepth());
            }
            else if (dest_file_type.equalsIgnoreCase("plink"))
            {
                ExportUtils.writeToPlink(gt_taxa_sites_filtered, dest_file, ' ');
            }
            else if (dest_file_type.equalsIgnoreCase("json"))
            {
                return writeToJSON(gt_taxa_sites_filtered, dest_file,splitAfterSite);
            }
            else
            {
                //System.err.println ("File format '" + dest_file_type + "' not recognized!");
                return "-1 -1\n" + "Output file format '" + dest_file_type + "' not recognized!";
            }
        }
        
        return return_Values;
    }
    private static String slice(String[] args)
    {
          //get parameters
          ArgsEngine  myArgsEngine = new ArgsEngine();
          myArgsEngine.add("-sf", "--sourcefile", true);
          myArgsEngine.add("-st", "--sourcefile-type", true);
          myArgsEngine.add("-df", "--destinationfile", true);
          myArgsEngine.add("-dt", "--destinationfile-type", true);
          myArgsEngine.add("-tf", "--taxa-file", true);
          myArgsEngine.add("-tl", "--taxa-list", true);
          myArgsEngine.add("-ch", "--chromosome", true);
          myArgsEngine.add("-start", "--chr-start", true);
          myArgsEngine.add("-end", "--chr-end", true);
          myArgsEngine.add("-bf", "--build-file", false);
          myArgsEngine.parse(args);
       
         boolean buildfile = false;
         String source_file = null;
         String source_file_type = null;
         String dest_file = null;
         String dest_file_type=null;
         String TaxaListFile=null;
         ArrayList<String> TaxaArrayList = new ArrayList<String>();
         String ChromosomeStr=null;
         int StartPhysicalPosition = 0;
         int EndPhysicalPosition = 0;

         if (myArgsEngine.getBoolean("-bf")) 
         {
              buildfile = true;      
         }
         
        if (myArgsEngine.getBoolean("-sf")) {
            source_file = myArgsEngine.getString("-sf");
        } else {
            printUsage("slice");
            throw new IllegalArgumentException("Please specify a source file (option -sf).");
        }
        
        if (myArgsEngine.getBoolean("-st")) {
            source_file_type = myArgsEngine.getString("-st");
        } else {
            printUsage("slice");
            throw new IllegalArgumentException("Please specify a source file type (option -st).");
        }
        
        if (buildfile)
        {
            if (myArgsEngine.getBoolean("-df")) {
                dest_file = myArgsEngine.getString("-df");
            } else {
                printUsage("slice");
                throw new IllegalArgumentException("Please specify a destination file (option -df).");
            }

            if (myArgsEngine.getBoolean("-dt")) {
                dest_file_type = myArgsEngine.getString("-dt");
            } else {
                printUsage("slice");
                throw new IllegalArgumentException("Please specify a destination file type (option -dt).");
            }
        }
        
        if (myArgsEngine.getBoolean("-tf")) 
        {
            TaxaListFile = myArgsEngine.getString("-tf");
            File outDirectory = new File(TaxaListFile);
            if (!outDirectory.isFile()) {
                printUsage("slice");
                throw new IllegalArgumentException("The taxa file you supplied (option -tf) is not a file: " + TaxaListFile);
            }
            //verify and create sub-taxalist
            //create taxa filtered genotype table
           
           try {
                BufferedReader br = new BufferedReader(new FileReader(TaxaListFile), 65536);

                String temp;
                int currLine = 0;
                while (((temp = br.readLine()) != null)) {
                    if (!temp.trim().isEmpty())
                    {
                        TaxaArrayList.add(temp.trim());
                        currLine++;
                    }

                }
            } catch (Exception e) {
                System.out.println("Couldn't open taxa file to read taxa list: " + e);
            }
        } 
        else if (myArgsEngine.getBoolean("-tl"))
        {
            String TaxaListStr = myArgsEngine.getString("-tf");
            if (TaxaListStr.equalsIgnoreCase("all"))
            {
                TaxaArrayList.add("all");
            }
            else
            {
                String[] wordList = TaxaListStr.split(";");
                TaxaArrayList.addAll(Arrays.asList(wordList));  
            }
        }
        else {
            printUsage("slice");
            throw new IllegalArgumentException("Please specify the file with taxa list (option -tf).");
        }
                               
        if (myArgsEngine.getBoolean("-ch")) {
            ChromosomeStr = myArgsEngine.getString("-ch");
        } else {
            if ((myArgsEngine.getBoolean("-start")) || (myArgsEngine.getBoolean("-end")))
            {
                printUsage("slice");
                throw new IllegalArgumentException("Please specify a chromosome name (option -ch).");
            }
            else
            {
                ChromosomeStr = "all";
            }
        }
        
        if (myArgsEngine.getBoolean("-start")) {
            StartPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-start"));
        } 
        
        if (myArgsEngine.getBoolean("-end")) {
            EndPhysicalPosition = Integer.parseInt(myArgsEngine.getString("-end"));
        }

         //create a genotypetable from source
        GenotypeTable source_gt_table = null;
        if (source_file_type.equalsIgnoreCase("hdf5"))
        {
            source_gt_table = GenotypeTableBuilder.getInstance(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("hapmap"))
        {
            source_gt_table = ImportUtils.readFromHapmap(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("vcf"))
        {
            source_gt_table = ImportUtils.readFromVCF(source_file, null, false);
        }

        

        
        GenotypeTable gt_taxa_filtered = null ; 
        
                
        
        //if TaxaListString=all or not specified, all taxa included
        if ((TaxaArrayList.get(0).equalsIgnoreCase("all")) || (TaxaArrayList.size()==0))
        {
            gt_taxa_filtered = source_gt_table;
        }
        //if TaxaListString specified, create subset
        else
        {
            TaxaList source_taxa_list = source_gt_table.taxa();
            TaxaList dest_taxa_list=null;
            TaxaListBuilder dest_taxa_list_builder = new TaxaListBuilder();

            ArrayList<String> unknow_taxa = new ArrayList<String>();
            for (int i =0; i<TaxaArrayList.size(); i++)
            {
               if (TaxaArrayList.get(i).matches(".*\\w.*"))
               {
                    int itaxa_index = source_taxa_list.indexOf(TaxaArrayList.get(i));
                    if (itaxa_index>=0)
                    {
                        dest_taxa_list_builder.add(source_taxa_list.get(itaxa_index));
                    }
                    else
                    {
                        unknow_taxa.add(TaxaArrayList.get(i));
                    }
               }
            }
            dest_taxa_list = dest_taxa_list_builder.build();
            if (unknow_taxa.size()>0)
            {
                String errmsg = "-1 -1\nUnknown taxa names:\n";
                for (int i=0; i<unknow_taxa.size(); i++)
                {
                    errmsg += unknow_taxa.get(i) + "\n";
                }
                return errmsg;
            }
            gt_taxa_filtered = FilterGenotypeTable.getInstance(source_gt_table, dest_taxa_list);
        }
        
        
        //create sites filtered genotype table
        GenotypeTable gt_taxa_sites_filtered = null;
        
        //create position filtered genotype table
        if (ChromosomeStr.toLowerCase().equalsIgnoreCase("all"))
        {
            gt_taxa_sites_filtered = gt_taxa_filtered;
        }
        else
        {
            Chromosome ch = source_gt_table.chromosome(ChromosomeStr);
            if (ch == null)
            {
                return "-1 -1\n" + "The chromosome name '" + ChromosomeStr +  "' is not recognized\n";
            }
            if ((StartPhysicalPosition==0) || (EndPhysicalPosition==0))
            {
                gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch);
            }
            else
            {
                gt_taxa_sites_filtered = FilterGenotypeTable.getInstance(gt_taxa_filtered, ch, StartPhysicalPosition, EndPhysicalPosition);
            }
            
        }
        String return_Values = "";
        if (gt_taxa_sites_filtered == null) 
        {
            return "-1 -1\nNo data after filtering!\n";
        }
        else
        {
            return_Values= gt_taxa_sites_filtered.taxa().size() + " " + gt_taxa_sites_filtered.positions().size();
        }
        if (buildfile)
        {
            if (dest_file_type.equalsIgnoreCase("hdf5"))
            {
                ExportUtils.writeGenotypeHDF5(gt_taxa_sites_filtered, dest_file, source_gt_table.hasDepth());
            }
            else if (dest_file_type.equalsIgnoreCase("hapmap"))
            {
                ExportUtils.writeToHapmap(gt_taxa_sites_filtered, dest_file);
            }
            else if (dest_file_type.equalsIgnoreCase("vcf"))
            {
                ExportUtils.writeToVCF(gt_taxa_sites_filtered, dest_file, source_gt_table.hasDepth());
            }
            else if (dest_file_type.equalsIgnoreCase("plink"))
            {
                ExportUtils.writeToPlink(gt_taxa_sites_filtered, dest_file, ' ');
            }
            else if (dest_file_type.equalsIgnoreCase("json"))
            {
                writeToJSON(gt_taxa_sites_filtered, dest_file,0);
            }
            else
            {
                //System.err.println ("File format '" + dest_file_type + "' not recognized!");
                return "-1 -1\n" + "Output file format '" + dest_file_type + "' not recognized!";
            }
        }
        
        return return_Values;
    }
    
    
    private static String[] get_taxa_list(String[] args)
    {
        ArgsEngine  myArgsEngine = new ArgsEngine();
          myArgsEngine.add("-sf", "--sourcefile", true);
          myArgsEngine.add("-st", "--sourcefile-type", true);
          myArgsEngine.parse(args);
       
         String source_file = null;
         String source_file_type = null;
         
        if (myArgsEngine.getBoolean("-sf")) {
            source_file = myArgsEngine.getString("-sf");
        } else {
            printUsage("get_taxa_list");
            throw new IllegalArgumentException("Please specify a source file (option -sf).");
        }
        
        if (myArgsEngine.getBoolean("-st")) {
            source_file_type = myArgsEngine.getString("-st");
        } else {
            printUsage("get_taxa_list");
            throw new IllegalArgumentException("Please specify a source file type (option -st).");
        }

        //create a genotypetable from source
        GenotypeTable source_gt_table = null;
        if (source_file_type.equalsIgnoreCase("hdf5"))
        {
            source_gt_table = GenotypeTableBuilder.getInstance(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("hapmap"))
        {
            source_gt_table = ImportUtils.readFromHapmap(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("vcf"))
        {
            source_gt_table = ImportUtils.readFromVCF(source_file, null, false);
        }
        
        TaxaList taxalist =source_gt_table.taxa();
        int taxaCount = taxalist.size();
        String[] taxa_array = new String[taxaCount];
        for (int i=0; i<taxaCount; i++)
        {
            taxa_array[i] = taxalist.taxaName(i);
        }
        return taxa_array;
    }
    
    private static String[] chr_marker_info(String[] args)
    {
        ArgsEngine  myArgsEngine = new ArgsEngine();
          myArgsEngine.add("-sf", "--sourcefile", true);
          myArgsEngine.add("-st", "--sourcefile-type", true);
          myArgsEngine.parse(args);
       
         String source_file = null;
         String source_file_type = null;
         
        if (myArgsEngine.getBoolean("-sf")) {
            source_file = myArgsEngine.getString("-sf");
        } else {
            printUsage("chr_marker_info");
            throw new IllegalArgumentException("Please specify a source file (option -sf).");
        }
        
        if (myArgsEngine.getBoolean("-st")) {
            source_file_type = myArgsEngine.getString("-st");
        } else {
            printUsage("chr_marker_info");
            throw new IllegalArgumentException("Please specify a source file type (option -st).");
        }

                 //create a genotypetable from source
        GenotypeTable source_gt_table = null;
        if (source_file_type.equalsIgnoreCase("hdf5"))
        {
            source_gt_table = GenotypeTableBuilder.getInstance(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("hapmap"))
        {
            source_gt_table = ImportUtils.readFromHapmap(source_file);
        }
        else if (source_file_type.equalsIgnoreCase("vcf"))
        {
            source_gt_table = ImportUtils.readFromVCF(source_file, null, false);
        }
        
        Chromosome[] chrlist =source_gt_table.chromosomes();
        int chrCount = chrlist.length;
        String[] chr_array = new String[chrCount];
        for (int i=0; i<chrCount; i++)
        {
            int[] firstlastpos = source_gt_table.firstLastSiteOfChromosome(chrlist[i]);
            int firstposition = source_gt_table.chromosomalPosition(firstlastpos[0]);
            int lastposition = source_gt_table.chromosomalPosition(firstlastpos[1]);
            chr_array[i] = chrlist[i].getName() + "\t" + source_gt_table.chromosomeSiteCount(chrlist[i]) +  "\t" + firstposition + "\t" + lastposition;
        }
        return chr_array;
    }
    
    // Only used by writeToJSON
    public static String appendIndexToPath(String fullFileName, int i){
    	String[] splitPath = fullFileName.split("/");
        splitPath[splitPath.length-1] = Integer.toString(i) + "_" + splitPath[splitPath.length-1];
        return String.join("/", splitPath);
    }
    
	/**
     * Write a GenotypeTable to JSON format
     * 
     * @param alignment genotype table
     * @param filename outfile name (will add ".json" if needed)
     * @param splitFileOnSite Splitting offset to split file after every n alleles
     * @return name of the final outfile with the appropriate suffix
     */
    public static String writeToJSON(GenotypeTable alignment, String filename, int splitFileOnSite){
    	int fileCount = 1;
    	
    	List<HashMap<String, String>> json_file = new ArrayList<HashMap<String, String>>();
    	ObjectMapper mapper = JsonFactory.create();
        String fullFileName = Utils.addSuffixIfNeeded(filename, ".json", new String[]{".json", ".JSON"});
         try {
             int numTaxa = alignment.numberOfTaxa();
             int numSites = alignment.numberOfSites();
         	 String[] results = new String[numTaxa];
             HashMap<String, String> hm = new HashMap<String, String>();
             for (int site = 0; site < numSites; site++) {
                 // STORE: SNP ID
            	 hm.put("rs#", alignment.siteName(site));
                 int[][] sortedAlleles = alignment.allelesSortedByFrequency(site); // which alleles are actually present among the genotypes
                 int numAlleles = sortedAlleles[0].length;
                 if (numAlleles == 0) {
                     //STORE: 'NA' if allele does not exist
                	 hm.put("allele", "NA");
                 } else if (numAlleles == 1) {
                	//STORE allele
                	 hm.put("allele",alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]));
                 } else { // if multiple
                	 String genotypeAtSite = alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]);
                     for (int allele = 1; allele < sortedAlleles[0].length; allele++) {
                         if (sortedAlleles[0][allele] != GenotypeTable.UNKNOWN_ALLELE) {
                             genotypeAtSite += "/";  
                             genotypeAtSite += alignment.genotypeAsString(site, (byte) sortedAlleles[0][allele]);
                             // STORE: multiple alleles
                             hm.put("allele",genotypeAtSite);
                         }
                     }
                 }
                //STORE: chrom_name, chrom_pos
                 hm.put("chrom_name", alignment.chromosomeName(site));
                 hm.put("chrom_pos", String.valueOf(alignment.chromosomalPosition(site)));
                 for (int taxa = 0; taxa < numTaxa; taxa++) {
                         String baseIUPAC = null;
                         try {
                             baseIUPAC = alignment.genotypeAsString(taxa, site);
                         } catch (Exception e) {
                             String[] b = alignment.genotypeAsStringArray(taxa, site);
                             myLogger.debug(e.getMessage(), e);
                             throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                         }
                         if ((baseIUPAC == null) || baseIUPAC.equals("?")) {
                             String[] b = alignment.genotypeAsStringArray(taxa, site);
                             throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                         }
                         results[taxa] = baseIUPAC;
                 }
                 hm.put("results", Arrays.toString(results));
                 json_file.add((HashMap<String, String>) hm.clone());
                 // Clear
                 hm.clear();
                 //split & write file
                 if (splitFileOnSite > 0 && site % splitFileOnSite == 0 && site > 0){
                	 mapper.writeValue(new FileOutputStream(appendIndexToPath(fullFileName,fileCount)), json_file);
                	 json_file.clear();
                	 fileCount++;
                 }
             }
             return Integer.toString(fileCount);
         } catch (Exception e) {
             myLogger.debug(e.getMessage(), e);
             throw new IllegalArgumentException("Error writing JSON file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
         } finally {
        	 try {
            	 System.out.println(Boon.toJson(json_file));
				mapper.writeValue(new FileOutputStream(appendIndexToPath(fullFileName,fileCount)), json_file);
			} catch (FileNotFoundException e) {
				 myLogger.debug(e.getMessage(), e);
	             throw new IllegalArgumentException("Error writing JSON file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
			}
         }
    }
    /**
     * 
     * @param menuName declare name of method to print usage
     * Prints out instructions on how to use the different methods
     * TODO: sliceNuke and sliceRangeNuke is missing.
     */
    private static void printUsage(String menuName)
    {
        if (menuName.equals("main"))
        {
            System.out.print("Usage:\tjava -jar Tassel_gt_server.jar [options]\n\n" 
                    +"Command:\tslice\tcreate a slice from the genotype file\n"
                    +"\t\tdbinfo\twrite the database information into files\n");
        }
        else if (menuName.equals("slice"))
        {
            System.out.print("Usage:\tjava -jar Tassel_gt_server.jar slice [options]\n\n" 
                    +"Options\t-bf\tbuild the file. If skip, only give the dimensions\n"
                    +"\t\t-sf\tsource file name\n"
                    +"\t\t-st\tsource file type (hdf5, vcf, or hapmap)\n"
                    +"\t\t-df\tdestination file name\n"
                    +"\t\t-dt\tdestination file type (hdf5, vcf, hapmap or plink)\n"
                    +"\t\t-tf\ttaxa file name\n"
                    +"\t\t-ch\tchromosome name\n"
                    +"\t\t-start\tphysical start position on the chromosome\n"
                   +"\t\t-end\tphysical end position on the chromosome\n"
                    );
        }
        else if (menuName.equals("get_taxa_list"))
        {
            System.out.print("Usage:\tjava -jar Tassel_gt_server.jar get_taxa_list [options]\n\n" 
                    +"Options\t-bf\tbuild the file. If skip, only give the dimensions\n"
                    +"\t\t-sf\tsource file name\n"
                    +"\t\t-st\tsource file type (hdf5, vcf, or hapmap)\n"
                    );
        }
        else if (menuName.equals("dbinfo"))
        {
            System.out.print("Usage:\tjava -jar Tassel_gt_server.jar dbinfo [options]\n\n" 
                    +"Options\t-bf\tbuild the file. If skip, only give the dimensions\n"
                    +"\t\t-sf\tsource file name\n"
                    +"\t\t-st\tsource file type (hdf5, vcf, or hapmap)\n"
                    +"\t\t-dbname\tdatabase name in one word\n"
                    +"\t\t-o\toutput directory\n"
                    );
        }
        else if (menuName.equals("chr_marker_info"))
        {
            System.out.print("Usage:\tjava -jar Tassel_gt_server.jar chr_marker_info [options]\n\n" 
                    +"Options\t-bf\tbuild the file. If skip, only give the dimensions\n"
                    +"\t\t-sf\tsource file name\n"
                    +"\t\t-st\tsource file type (hdf5, vcf, or hapmap)\n"
                    );
        }
        else
        {
            printUsage("main");
        }
    }
        
}

