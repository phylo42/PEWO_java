
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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import jplace.JplacerLoader;
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
                String readLabel=currentJPlaceFile.getFileName().toString().split("\\.")[1]; //SECOND element in filename in pplacer ex: RAxML_portableTree.R0_nx110_la_r150.jplace
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
            
            System.out.println("##############");
            System.out.println("## PPL");
            File PPLxDir=new File(dg.workDir+File.separator+"PPLx");
            //load EPA jplace results
            ArrayList<File> PPLResults=new ArrayList<>();
            List<Path> PPLJPlaceFiles = Files.find(PPLxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

            //for each jplace file, calculate node dist to expected placement
            for (int i = 0; i < PPLJPlaceFiles.size(); i++) {
                Path currentJPlaceFile = PPLJPlaceFiles.get(i);
                String experimentLabel=currentJPlaceFile.getParent().getFileName().toString();
                String readLabel=currentJPlaceFile.getFileName().toString().split("\\.")[0]; //FIRST element in filename in pplacer ex: R0_nx110_la_r150.aln.jplace
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
                    bw.append("PPL;"+experimentLabel+";;;"+readLabel+";"+name+";"+nodeDistance+"\n");
                }
            
            }
            
            System.out.println("##############");
            System.out.println("## RAP");
            File DxDir=new File(dg.workDir+File.separator+"Dx");
            //load EPA jplace results
            ArrayList<File> RAPResults=new ArrayList<>();
            List<Path> RAPPJPlaceFiles = Files.find(DxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

            //for each jplace file, calculate node dist to expected placement
            for (int i = 0; i < RAPPJPlaceFiles.size(); i++) {
                Path currentJPlaceFile = RAPPJPlaceFiles.get(i);
                System.out.println(currentJPlaceFile.toAbsolutePath());
                //k and alpha
                String kAlphaLabel=currentJPlaceFile.getParent().getParent().getFileName().toString(); //2 times get parent  kx_ax/logs/jplace
                String[] data =kAlphaLabel.split("_");
                int k=Integer.parseInt(data[0].substring(1));
                float alpha=Float.parseFloat(data[1].substring(1));
                //experiment
                String experimentLabel=currentJPlaceFile.getParent().getParent().getParent().getFileName().toString(); //3 times get parent  Ax_nxx_xxx/kx_ax/logs/jplace
                String readLabel=currentJPlaceFile.getFileName().toString().split("\\.")[0]; //FIRST element in filename in pplacer ex: R0_nx110_la_r150.aln.jplace
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
                    bw.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+readLabel+";"+name+";"+nodeDistance+"\n");
                }
            
            } 
            
            bw.close();
            
            System.exit(0);
            
            
        
            
            
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
