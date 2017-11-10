
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import jplace.JplacerLoader;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class CompareLikelihoodRatios {

    public static void main(String[] args) throws IOException {
        
        File RAPFile=new File("/home/ben/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML/logs/placements_mod_p4z1r36_query_only2.fasta_union.jplace");
        File EPAFile=new File("/home/ben/Dropbox/viromeplacer/test_datasets/pplacer_epa_big_dataset/epa/RAxML_portableTree.EPA.jplace");
        File PPLFile=new File("/media/ben/STOCK/SOFTWARE/pplacer-Linux-v1.1.alpha18-2-gcb55169/fhcrc-microbiome-demo-730d268/p4z1r36.jplace");

        JplacerLoader jlRAP=new JplacerLoader(RAPFile);
        JplacerLoader jlEPA=new JplacerLoader(EPAFile);
        JplacerLoader jlPPL=new JplacerLoader(PPLFile);
        
        File comparisonFile=new File("/media/ben/STOCK/DATA/viromeplacer/accu_tests/R_analysis/comparison_weightratio.csv");
        BufferedWriter writer = Files.newBufferedWriter(comparisonFile.toPath(), StandardCharsets.UTF_8);   
        writer.append("read\tsoftware\twr\n");
        //writer.append("read\tidx\tRAPwr\tEPAwr\tPPLwr\n");
        
        
        NumberFormat nf = NumberFormat.getNumberInstance();
        nf.setMaximumFractionDigits(12);
        
        System.out.println(jlRAP.getWeightRatios().get("GLKT0ZE01BPHOX"));
        
        
        int readCounter=0;
        for (Iterator<String> seqIterator = jlRAP.getWeightRatios().keySet().iterator(); seqIterator.hasNext();) {
            String seqName = seqIterator.next();
            //System.out.println("seqName:"+seqName);
            HashMap<String, ArrayList<Double>> RAPWeightRatios = jlRAP.getWeightRatios();
            HashMap<String, ArrayList<Double>> EPAWeightRatios = jlEPA.getWeightRatios();
            HashMap<String, ArrayList<Double>> PPLWeightRatios = jlPPL.getWeightRatios();
            
//            int placement=0;
//            
//            int max=RAPWeightRatios.get(seqName).size();
//            if (EPAWeightRatios.get(seqName).size()>max)
//                max=EPAWeightRatios.get(seqName).size();
//            if (PPLWeightRatios.get(seqName).size()>max)
//                max=PPLWeightRatios.get(seqName).size();
//            
//            for (int i = 0; i < 1; i++) {
//                writer.append(readCounter+"\t"+i+"\t");
//                if (i<RAPWeightRatios.get(seqName).size())
//                    writer.append(nf.format(RAPWeightRatios.get(seqName).get(i))+"\t");
//                else 
//                    writer.append("\t");
//                if (i<EPAWeightRatios.get(seqName).size())
//                    writer.append(nf.format(EPAWeightRatios.get(seqName).get(i))+"\t");
//                else 
//                    writer.append("\t");
//                if (i<PPLWeightRatios.get(seqName).size())
//                    writer.append(nf.format(PPLWeightRatios.get(seqName).get(i))+"\n");
//                else 
//                    writer.append("\n");
//            }



            for (int i = 0; i < RAPWeightRatios.get(seqName).size(); i++) {
                writer.append(readCounter+"\tRAP\t"+nf.format(RAPWeightRatios.get(seqName).get(i))+"\n");
            }
            for (int i = 0; i < EPAWeightRatios.get(seqName).size(); i++) {
                writer.append(readCounter+"\tEPA\t"+nf.format(EPAWeightRatios.get(seqName).get(i))+"\n");
            }
            for (int i = 0; i < PPLWeightRatios.get(seqName).size(); i++) {
                writer.append(readCounter+"\tPPL\t"+nf.format(PPLWeightRatios.get(seqName).get(i))+"\n");
            }
            
            readCounter++;
        }
        
        writer.close();

        System.out.println("DONE");
        
    }
    
    
    

    
    
    
    
}
