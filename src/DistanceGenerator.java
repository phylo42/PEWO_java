
import alignement.Alignment;
import core.States;
import dtx.Dtx;
import etc.Infos;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
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

    
    public static void main(String[] args) {
        
        ObjectInputStream ois=null;
        try {
            
            System.out.println("ARGS: workDir");
            
            //launch
            DistanceGenerator dg=new DistanceGenerator();
            
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                dg.workDir=new File(args[0]);
                System.out.println("workDir: "+dg.workDir);
            }  
            
            //expected placement
            File expPLaceFile=new File(dg.workDir+File.separator+"expected_placements.bin");
            //Dtx
            File DtxFile=new File(dg.workDir+File.separator+"Dtx.csv");
           

            //load the expected placement
            ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(expPLaceFile),4096));
            Infos.println("Loading NxIndex");
            HashMap<Integer,Integer> NxIndex = (HashMap<Integer,Integer>)ois.readObject();
            Infos.println("Loading pruned trees");
            ArrayList<ArrayList<Integer>> expectedPlacements = (ArrayList<ArrayList<Integer>>)ois.readObject();
            Infos.println("Loading expected placements");
            ArrayList<PhyloTree> prunedTrees = (ArrayList<PhyloTree>)ois.readObject();
            
            //load Dtx
            Dtx Dtx=new Dtx(DtxFile);
            System.out.println(Dtx);
            
            //prepare a nice CSV file in which all data will be saved
            File csvResult=new File(dg.workDir+File.separator+"results.csv");
            BufferedWriter bw=new BufferedWriter(new FileWriter(csvResult));
            //header
            bw.append("software;Ax;k;alpha;Rx;read;node_dist\n");
            
            
            
            System.out.println("##############");
            System.out.println("## EPA");
            File EPAxDir=new File(dg.workDir+File.separator+"EPAx");
            //load EPA jplace results
            ArrayList<File> EPAResults=new ArrayList<>();
            List<Path> EPAJPlaceFiles = Files.find(EPAxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

            //for each jplace file, calculate node dist to expected placement
            for (int i = 0; i < EPAJPlaceFiles.size(); i++) {
                Path currentJPlaceFile = EPAJPlaceFiles.get(i);
                String experimentLabel=currentJPlaceFile.getParent().getFileName().toString();
                String readLabel=currentJPlaceFile.getFileName().toString().split("\\.")[1];
                int pruningNumber=Integer.parseInt(experimentLabel.split("_")[0].substring(1));
                int prunedNodeId=Integer.parseInt(experimentLabel.split("_")[1].substring(2)); //i.e. Nx
                int expectedPlacementIndex=NxIndex.get(prunedNodeId);
                System.out.println("experimentLabel:"+experimentLabel+" read:"+readLabel);
                
                //load tree and expectedPlacement related to this Ax
                PhyloTree experimentTree=prunedTrees.get(expectedPlacementIndex);
                ArrayList<Integer> experimentPlacements=expectedPlacements.get(expectedPlacementIndex);
                //System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
                //System.out.println("experimentTree best placement(s):"+experimentPlacements);
                
                JplacerLoader EPAJplace=new JplacerLoader(currentJPlaceFile.toFile());
                //map EPA jplace to experimentTree
                HashMap<Integer, Integer> mapEPANodes = EPAJplace.getTree().mapNodes(experimentTree);
                //retrieve best placements
                HashMap<String, Integer> EPABestPlacements = EPAJplace.getBestPlacements();


                //System.out.println("mapEPANodes:"+mapEPANodes);
                //System.out.println("EPABestPlacements:"+EPABestPlacements);

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

                    int nodeDistance = Dtx.getNodeDistance(prunedNodeId, experimentTreeNodeId);

                    System.out.println(name+" -> nodeDistance:"+nodeDistance);
                    //software;Ax;k;alpha;Rx;read;node_dist
                    bw.append("EPA;"+experimentLabel+";;;"+readLabel+";"+name+";"+nodeDistance+"\n");
                }
            
            }
                
            
            bw.close();
            
            System.exit(0);
            
            //TODO: taxit not working on server for now, need python libs and
            //environements !
            
//            System.out.println("##############");
//            System.out.println("## PPL");
//            
//            File test2=new File("/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set/PPLx/A0_nx0_laW/R0_nx0_laW_r3.aln.jplace");
//            JplacerLoader PPLJplace=new JplacerLoader(test2);
//            HashMap<Integer, Integer> mapPPLNodes = PPLJplace.getTree().mapNodes(experimentTree);
//            HashMap<String, Integer> PPLBestPlacements = PPLJplace.getBestPlacements();
//
//            System.out.println("mapPPLNodes:"+mapPPLNodes);
//            System.out.println("PPLBestPlacements:"+PPLBestPlacements);
//            for (Iterator<String> iterator = PPLBestPlacements.keySet().iterator(); iterator.hasNext();) {
//                String name = iterator.next();
//                //get best placement as the nodeId of the phylotree generated 
//                //during jplace parsing
//                Integer jplacePhyloTreeNodeId = PPLBestPlacements.get(name);
//                //get its equivalent nodeId in the phylotree loaded from the 
//                //expected_placements.bin
//                Integer experimentTreeNodeId = mapPPLNodes.get(jplacePhyloTreeNodeId);
//                //calculate the distance between these 2 nodeIds
//                //i.e. use the DTx and D'Tx matrices
//                
//                //TODO: add the Nx ids in the expected placement binary
//                
//                int nodeDistance = Dtx.getNodeDistance(0, experimentTreeNodeId);
//                
//                System.out.println(name+" -> nodeDistance:"+nodeDistance);
//                
//            }
            
            
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
