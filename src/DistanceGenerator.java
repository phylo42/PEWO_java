
import alignement.Alignment;
import core.States;
import dtx.Dtx;
import etc.Infos;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import jplace.JplacerLoader;
import main_v2.SessionNext_v2;
import tree.ExtendedTree;
import tree.PhyloTree;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class DistanceGenerator {
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");
    //expected placement
    File expPLaceFile=new File(workDir+File.separator+"expected_placements.bin");
    //Dtx
    File DtxFile=new File(workDir+File.separator+"Dtx.csv");
    
    public static void main(String[] args) {
        
        ObjectInputStream ois=null;
        try {
            DistanceGenerator dg=new DistanceGenerator();
            //first test on one case (A0), then change to use everything in Ax.
            
            //load EPA jplace result
            File test=new File("/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set/EPAx/A0_nx0_laW/RAxML_portableTree.R0_nx0_laW_r3.jplace");
            JplacerLoader EPAJplace=new JplacerLoader(test);
            File test2=new File("/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set/PPLx/A0_nx0_laW/R0_nx0_laW_r3.aln.jplace");
            JplacerLoader PPLJplace=new JplacerLoader(test2);
            //load the expected placement
            ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(dg.expPLaceFile),4096));
            Infos.println("Loading pruned trees");
            ArrayList<ArrayList<Integer>> expectedPlacements = (ArrayList<ArrayList<Integer>>)ois.readObject();
            Infos.println("Loading expected placements");
            ArrayList<PhyloTree> prunedTrees = (ArrayList<PhyloTree>)ois.readObject();
            
            PhyloTree experimentTree=prunedTrees.get(0);
            ArrayList<Integer> experimentPlacements=expectedPlacements.get(0);
            System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
            System.out.println("experimentTree best placement(s):"+experimentPlacements);
            
            //map EPA jplace to experimentTree
            HashMap<Integer, Integer> mapEPANodes = EPAJplace.getTree().mapNodes(experimentTree);
            HashMap<Integer, Integer> mapPPLNodes = PPLJplace.getTree().mapNodes(experimentTree);
            
            
            HashMap<String, Integer> EPABestPlacements = EPAJplace.getBestPlacements();
            HashMap<String, Integer> PPLBestPlacements = PPLJplace.getBestPlacements();
            
            //load Dtx
            Dtx Dtx=new Dtx(dg.DtxFile);
            System.out.println(Dtx);
            
            System.out.println("##############");
            System.out.println("## EPA");
            System.out.println("mapEPANodes:"+mapEPANodes);
            System.out.println("EPABestPlacements:"+EPABestPlacements);
            
            for (Iterator<String> iterator = EPABestPlacements.keySet().iterator(); iterator.hasNext();) {
                String name = iterator.next();
                //get best placement as the nodeId of the phylotree generated 
                //during jplace parsing
                Integer jplacePhyloTreeNodeId = EPABestPlacements.get(name);
                //get its equivalent nodeId in the phylotree loaded from the 
                //expected_placements.bin
                Integer experimentTreeNodeId = mapEPANodes.get(jplacePhyloTreeNodeId);
                //calculate the distance between these 2 nodeIds
                //i.e. use the DTx and D'Tx matrices
                
                //TODO: add the Nx ids in the expected placement binary
                
                int nodeDistance = Dtx.getNodeDistance(0, experimentTreeNodeId);
                
                System.out.println(name+" -> nodeDistance:"+nodeDistance);
                
            }
            
            System.out.println("##############");
            System.out.println("## PPL");
            System.out.println("mapPPLNodes:"+mapPPLNodes);
            System.out.println("PPLBestPlacements:"+PPLBestPlacements);
            for (Iterator<String> iterator = PPLBestPlacements.keySet().iterator(); iterator.hasNext();) {
                String name = iterator.next();
                //get best placement as the nodeId of the phylotree generated 
                //during jplace parsing
                Integer jplacePhyloTreeNodeId = PPLBestPlacements.get(name);
                //get its equivalent nodeId in the phylotree loaded from the 
                //expected_placements.bin
                Integer experimentTreeNodeId = mapPPLNodes.get(jplacePhyloTreeNodeId);
                //calculate the distance between these 2 nodeIds
                //i.e. use the DTx and D'Tx matrices
                
                //TODO: add the Nx ids in the expected placement binary
                
                int nodeDistance = Dtx.getNodeDistance(0, experimentTreeNodeId);
                
                System.out.println(name+" -> nodeDistance:"+nodeDistance);
                
            }
            
            
            //TODO: chack that these mappings are OK !
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DistanceGenerator.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException | ClassNotFoundException ex) {
            Logger.getLogger(DistanceGenerator.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                ois.close();
            } catch (IOException ex) {
                Logger.getLogger(DistanceGenerator.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        
        
        
    }
    
}
