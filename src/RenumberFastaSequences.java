
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class RenumberFastaSequences {
    
    public static void main(String[] args) {
        
        try {
            File fasta=new File("/media/ben/STOCK/DATA/viromeplacer/speed_tests/EMP_92_studies_10000000.fas_pplacer16S.aln.fa");
            
            BufferedReader br = Files.newBufferedReader(fasta.toPath());
            
            BufferedWriter bw = Files.newBufferedWriter(new File("/media/ben/STOCK/DATA/viromeplacer/speed_tests/alignment_10000000.fasta").toPath(),StandardCharsets.UTF_8);
            
            String line=null;
            int readCounter=0;
            while ((line=br.readLine())!=null) {
                if (readCounter%10000==0) {System.out.println(readCounter);}
                if (line.startsWith(">")) {
                    if (line.contains("_")) {
                        line=">r_"+readCounter;
                        readCounter++;
                    }
                }
                bw.append(line);
                bw.newLine();   
            }
            
            System.out.println("Done: "+readCounter);
            
            br.close();
            bw.close();
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(RenumberFastaSequences.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
    }
    
    
}
