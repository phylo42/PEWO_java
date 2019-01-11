
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.logging.Level;
import java.util.logging.Logger;
import jplace.JplacerLoader;
import tree.NewickWriter;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class RerootEPANGJplace {

    public static void main(String[] args) {
        
        try {
            File f=new File("/home/ben/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/GREEN85/epa_result_preplacement.jplace");
            File out=new File("/home/ben/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/GREEN85/test");
            BufferedWriter bw = Files.newBufferedWriter(out.toPath());
            
            JplacerLoader jpl=new JplacerLoader(f, true);
            
            NewickWriter nw=new NewickWriter(out);
            nw.writeNewickTree(jpl.getTree(), true, true, true, false);
            nw.close();
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(RerootEPANGJplace.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    
    
    
}
