
import dtx.Dtx;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import jplace.JplacerLoader;
import jplace.JplacerLoader.Placement;
import tree.PhyloNode;
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

    //which component to score ?
    boolean doEPA=true;
    boolean doPPL=true;
    boolean doRAP=true;
    boolean doEPANG=true;
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose","1");
        
        HashMap<String,Integer> PPLResults=new HashMap<String,Integer>(); //map(experimentlabel_read_name)=ARTree_nodeid
        HashMap<String,Integer> EPAResults=new HashMap<String,Integer>(); //map(experimentlabel_read_name)=node_dist
        ObjectInputStream ois=null;
        
        try {
            
            System.out.println("ARGS: workDir doEPA[1/0] doEPANG[0/1] doPPL[1/0] doRAP[0/1]");
            
            //launch
            DistanceGenerator dg=new DistanceGenerator();
            
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                dg.workDir=new File(args[0]);
                System.out.println("workDir: "+dg.workDir);
                if (Integer.parseInt(args[1])<1) { dg.doEPA=false; }
                if (Integer.parseInt(args[2])<1) { dg.doEPANG=false; }
                if (Integer.parseInt(args[3])<1) { dg.doPPL=false; }
                if (Integer.parseInt(args[4])<1) { dg.doRAP=false; }
                
            }  
            
            //expected placement
            File expPLaceFile=new File(dg.workDir+File.separator+"expected_placements.bin");
            //Dtx
            File DtxFile=new File(dg.workDir+File.separator+"Dtx.csv");
           

            //load the expected placement
            ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(expPLaceFile),4096));
            System.out.println("Loading NxIndex");
            HashMap<Integer,Integer> NxIndex = (HashMap<Integer,Integer>)ois.readObject();
            System.out.println("Loading expected placements");
            ArrayList<ArrayList<Integer>> expectedPlacementsNodeIds = (ArrayList<ArrayList<Integer>>)ois.readObject();
            System.out.println("Loading trees");
            ArrayList<PhyloTree> experimentTrees = (ArrayList<PhyloTree>)ois.readObject();
            
            System.out.println("################################################");
            System.out.println("NxIndex="+NxIndex);
            System.out.println("expectedPlacements="+expectedPlacementsNodeIds);
            System.out.println("prunedTrees="+experimentTrees.size());
            for (int i = 0; i < experimentTrees.size(); i++) {
                PhyloTree get = experimentTrees.get(i);
                System.out.println(i+"th tree size test:"+get.getNodeCount());
            }
            System.out.println("################################################");
            //load Dtx
            System.out.println("Loading Dtx matrix...");
            Dtx Dtx=new Dtx(DtxFile);
            //System.out.println(Dtx);
            
            //prepare a nice CSV file in which all data will be saved
            Path csvResult=Paths.get(dg.workDir.getAbsolutePath(),"results.csv");
            BufferedWriter bw= Files.newBufferedWriter(csvResult);
            //header
            bw.append("software;Ax;k;alpha;dbSize;Rx;read;readSize;node_dist;readStart;readEnd\n");
            
            //prepare a second CSV file confronting the methods
            Path csvResult2=Paths.get(dg.workDir.getAbsolutePath(),"results2.csv");
            BufferedWriter bw2= Files.newBufferedWriter(csvResult2);
            //header
            bw2.append("experiment;read;k;alpha;dbSize;readSize;ndist_EPA;ndist_PPL;ndist_RAP\n");
            
            //prepare a third CSV for P,B,E,L,R test
            Path csvResult3=Paths.get(dg.workDir.getAbsolutePath(),"results3.csv");
            BufferedWriter bw3= Files.newBufferedWriter(csvResult3);
            //header
            bw3.append("software;Ax;k;alpha;dbSize;Rx;read;readSize;category\n");
            
            
            if (dg.doEPA) {
                System.out.println("##############");
                System.out.println("## EPA");
                File EPAxDir=new File(dg.workDir+File.separator+"EPAx");
                //load EPA jplace results
                List<Path> EPAJPlaceFiles = Files.find(EPAxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

                //for each jplace file, calculate node dist to expected placement
                for (int i = 0; i < EPAJPlaceFiles.size(); i++) {
                    Path currentJPlaceFile = EPAJPlaceFiles.get(i);
                    String experimentLabel=currentJPlaceFile.getParent().getFileName().toString();
                    //SECOND element in filename in pplacer ex: RAxML_portableTree. R63_nx63_la1a.HCV_1_r200 .jplace
                    String readLabel=currentJPlaceFile.getFileName().toString().split("\\.jplace$")[0].substring(19); 
                    String[]elts=readLabel.split("_");
                    //ex: R63 nx63 la1a.HCV 1 r200
                    String readSize=elts[elts.length-1].substring(1);
                    int pruningNumber=Integer.parseInt(experimentLabel.split("_")[0].substring(1));
                    int prunedNodeId=Integer.parseInt(experimentLabel.split("_")[1].substring(2)); //i.e. Nx
                    int expectedPlacementIndex=NxIndex.get(prunedNodeId);
                    System.out.println("--------------------------------------");
                    System.out.println("experimentLabel:"+experimentLabel+" read:"+readLabel);

                    //load tree and expectedPlacement related to this Ax
                    PhyloTree experimentTree=experimentTrees.get(expectedPlacementIndex);
                    experimentTree.initIndexes();
                    ArrayList<Integer> experimentPlacements=expectedPlacementsNodeIds.get(expectedPlacementIndex);
                    //System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("experimentTree best placement(s):"+expectedPlacementNodeIds);

                    
                    JplacerLoader EPAJplace=new JplacerLoader(currentJPlaceFile.toFile(), false);
                    //System.out.println("EPANGJplace tree: "+EPANGJplace.getTree());
                    //System.out.println("experimentTree: "+experimentTree);
                    //System.out.println("RAPJplace tree nodes ids by DFS:"+RAPJplace.getTree().getNodeIdsByDFS());
                    //System.out.println("experimentTree nodes ids by DFS:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("RAPJplace tree nodes by DFS:"+RAPJplace.getTree().getNodeIdsByDFS().stream().map((id)->RAPJplace.getTree().getById(id)).peek((id)-> System.out.println(id)).count());
                    //System.out.println("experimentTree nodes by DFS:"+experimentTree.getNodeIdsByDFS().stream().map((id)->experimentTree.getById(id)).peek((id)-> System.out.println(id)).count());
                    //System.out.println("JPlace best placements:"+RAPJplace.getPlacements());
//                    System.out.println("POSTERIOR");
//                    testPosteriorDFS(experimentTree.getRoot());
//                    System.out.println("ANTERIOR");
//                    testAnteriorDFS(experimentTree.getRoot());
                    
                    
                    if (EPAJplace.getTree().getNodeCount()!=experimentTree.getNodeCount()) {
                        System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                        System.out.println("They do not include the same trees for the same Nx experiment.");
                        System.exit(1);
                    }
                    
                    
                    
                    
                    
                    //map EPA jplace to experimentTree
                    HashMap<Integer, Integer> mapEPANodes = EPAJplace.getTree().mapNodes(experimentTree);
                    //retrieve best placements
                    HashMap<String, ArrayList<Placement>> EPABestPlacements = EPAJplace.getPlacements();


                    //System.out.println("mapPPLNodes:"+mapPPLNodes);
                    //System.out.println("RAPBestPlacements:"+RAPBestPlacements);

                    for (Iterator<String> iterator = EPABestPlacements.keySet().iterator(); iterator.hasNext();) {
                        String name = iterator.next();
                        //get best placement as the nodeId of the phylotree generated 
                        //during jplace parsing
                        Integer jplacePhyloTreeNodeId = EPABestPlacements.get(name).get(0).getNodeId();
                        //get its equivalent nodeId in the phylotree loaded from the 
                        //expected_placements.bin
                        Integer experimentTreeNodeId = mapEPANodes.get(jplacePhyloTreeNodeId);
                        //calculate the distance between these 2 nodeIds
                        //i.e. use the DTx and D'Tx matrices
                        int nodeDistance = Dtx.getNodeDistance(prunedNodeId, experimentTreeNodeId);
                        
                        //got coordinates of placed read
                        String[] readInfos=name.split("_");
                        int readStart=Integer.decode(readInfos[readInfos.length-2]);
                        int readEnd=Integer.decode(readInfos[readInfos.length-1]);

                        //System.out.println(name+" -> nodeDistance:"+nodeDistance);
                        EPAResults.put(experimentLabel+":"+name, nodeDistance);
                        //software;Ax;k;alpha;Rx;read;node_dist
                        bw.append("EPA;"+experimentLabel+";;;;"+readLabel+";"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
                    }

                }
            }    
            
            if (dg.doEPANG) {
                System.out.println("##############");
                System.out.println("## EPA-ng");
                File EPANGxDir=new File(dg.workDir+File.separator+"EPANGx");
                //load EPA jplace results
                List<Path> EPANGJPlaceFiles = Files.find(EPANGxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

                //for each jplace file, calculate node dist to expected placement
                for (int i = 0; i < EPANGJPlaceFiles.size(); i++) {
                    Path currentJPlaceFile = EPANGJPlaceFiles.get(i);
                    String experimentLabel=currentJPlaceFile.getParent().getFileName().toString();
                    //FIRST  element in filename in epa-ng ex: R0_nx3_la95_r150. aln. fasta_queries_epa_result. jplace
                    String readLabel=currentJPlaceFile.getFileName().toString().split("\\.aln\\.fasta_queries_epa_result\\.jplace$")[0]; 
                    String[]elts=readLabel.split("_");
                    //ex: R63 nx63 la1a.HCV 1 r200
                    String readSize=elts[elts.length-1].substring(1);
                    int pruningNumber=Integer.parseInt(experimentLabel.split("_")[0].substring(1));
                    int prunedNodeId=Integer.parseInt(experimentLabel.split("_")[1].substring(2)); //i.e. Nx
                    int expectedPlacementIndex=NxIndex.get(prunedNodeId);
                    System.out.println("--------------------------------------");
                    System.out.println("experimentLabel:"+experimentLabel+" read:"+readLabel);

                    //load tree and expectedPlacement related to this Ax
                    PhyloTree experimentTree=experimentTrees.get(expectedPlacementIndex);
                    experimentTree.initIndexes();
                    ArrayList<Integer> experimentPlacements=expectedPlacementsNodeIds.get(expectedPlacementIndex);
                    //System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("experimentTree best placement(s):"+expectedPlacementNodeIds);

                    
                    JplacerLoader EPANGJplace=new JplacerLoader(currentJPlaceFile.toFile(), true);
                    //System.out.println("EPANGJplace tree: "+EPANGJplace.getTree());
                    //System.out.println("experimentTree: "+experimentTree);
                    //System.out.println("RAPJplace tree nodes ids by DFS:"+RAPJplace.getTree().getNodeIdsByDFS());
                    //System.out.println("experimentTree nodes ids by DFS:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("RAPJplace tree nodes by DFS:"+RAPJplace.getTree().getNodeIdsByDFS().stream().map((id)->RAPJplace.getTree().getById(id)).peek((id)-> System.out.println(id)).count());
                    //System.out.println("experimentTree nodes by DFS:"+experimentTree.getNodeIdsByDFS().stream().map((id)->experimentTree.getById(id)).peek((id)-> System.out.println(id)).count());
                    //System.out.println("JPlace best placements:"+RAPJplace.getPlacements());
//                    System.out.println("POSTERIOR");
//                    testPosteriorDFS(experimentTree.getRoot());
//                    System.out.println("ANTERIOR");
//                    testAnteriorDFS(experimentTree.getRoot());

                    //current version of EPA-ng unroots the input tree
                    //this needs to be corrected
                    
                    
                    
                    if (EPANGJplace.getTree().getNodeCount()!=experimentTree.getNodeCount()) {
                        System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                        System.out.println("They do not include the same trees for the same Nx experiment.");
                        System.exit(1);
                    }
                    
                    //map EPA jplace to experimentTree
                    HashMap<Integer, Integer> mapEPANGNodes = EPANGJplace.getTree().mapNodes(experimentTree);
                    //retrieve best placements
                    HashMap<String, ArrayList<Placement>> EPANGBestPlacements = EPANGJplace.getPlacements();


                    //System.out.println("mapPPLNodes:"+mapPPLNodes);
                    //System.out.println("RAPBestPlacements:"+RAPBestPlacements);

                    for (Iterator<String> iterator = EPANGBestPlacements.keySet().iterator(); iterator.hasNext();) {
                        String name = iterator.next();
                        //get best placement as the nodeId of the phylotree generated 
                        //during jplace parsing
                        Integer jplacePhyloTreeNodeId = EPANGBestPlacements.get(name).get(0).getNodeId();
                        //get its equivalent nodeId in the phylotree loaded from the 
                        //expected_placements.bin
                        Integer experimentTreeNodeId = mapEPANGNodes.get(jplacePhyloTreeNodeId);
                        //calculate the distance between these 2 nodeIds
                        //i.e. use the DTx and D'Tx matrices
                        int nodeDistance = Dtx.getNodeDistance(prunedNodeId, experimentTreeNodeId);
                        
                        //got coordinates of placed read
                        String[] readInfos=name.split("_");
                        int readStart=Integer.decode(readInfos[readInfos.length-2]);
                        int readEnd=Integer.decode(readInfos[readInfos.length-1]);

                        //System.out.println(name+" -> nodeDistance:"+nodeDistance);
                        EPAResults.put(experimentLabel+":"+name, nodeDistance);
                        //software;Ax;k;alpha;Rx;read;node_dist
                        bw.append("EPA-ng;"+experimentLabel+";;;;"+readLabel+";"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
                    }

                }
            }    
            
            
            if (dg.doPPL) {
                
                System.out.println("##############");
                System.out.println("## PPL");
                File PPLxDir=new File(dg.workDir+File.separator+"PPLx");
                //load EPA jplace results
                List<Path> PPLJPlaceFiles = Files.find(PPLxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

                //for each jplace file, calculate node dist to expected placement
                for (int i = 0; i < PPLJPlaceFiles.size(); i++) {
                    Path currentJPlaceFile = PPLJPlaceFiles.get(i);
                    String experimentLabel=currentJPlaceFile.getParent().getFileName().toString();
                    String readLabel=currentJPlaceFile.getFileName().toString().split("\\.aln.jplace")[0]; //FIRST element in filename in pplacer ex: R0_nx110_la_r150.aln.jplace
                    String[]elts=readLabel.split("_");
                    String readSize=elts[elts.length-1].substring(1);
                    int pruningNumber=Integer.parseInt(experimentLabel.split("_")[0].substring(1));
                    int prunedNodeId=Integer.parseInt(experimentLabel.split("_")[1].substring(2)); //i.e. Nx
                    int expectedPlacementIndex=NxIndex.get(prunedNodeId);
                    System.out.println("experimentLabel:"+experimentLabel+" read:"+readLabel);

                    //load tree and expectedPlacement related to this Ax
                    PhyloTree experimentTree=experimentTrees.get(expectedPlacementIndex);
                    ArrayList<Integer> expectedPlacementNodeIds=expectedPlacementsNodeIds.get(expectedPlacementIndex);
                    //System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("experimentTree best placement(s):"+expectedPlacementNodeIds);

                    JplacerLoader PPLJplace=new JplacerLoader(currentJPlaceFile.toFile(), false);
                    //System.out.println("PPLJplace tree: "+PPLJplace.getTree());
                    //System.out.println("experimentTree: "+experimentTree);
                    if (PPLJplace.getTree().getNodeCount()!=experimentTree.getNodeCount()) {
                        System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                        System.out.println("They do not include the same trees for the same Nx experiment.");
                        System.exit(1);
                    }
                    
                    //map EPA jplace to experimentTree
                    HashMap<Integer, Integer> mapPPLNodes = PPLJplace.getTree().mapNodes(experimentTree);
                    //retrieve best placements
                    HashMap<String, ArrayList<Placement>> EPABestPlacements = PPLJplace.getPlacements();


                    //System.out.println("mapPPLNodes:"+mapPPLNodes);
                    //System.out.println("RAPBestPlacements:"+RAPBestPlacements);

                    for (Iterator<String> iterator = EPABestPlacements.keySet().iterator(); iterator.hasNext();) {
                        //query itself
                        String name = iterator.next();
                        //get best placement as the nodeId of the phylotree generated 
                        //during jplace parsing
                        Integer jplacePhyloTreeNodeId = EPABestPlacements.get(name).get(0).getNodeId();
                        //get its equivalent nodeId in the phylotree loaded from the 
                        //expected_placements.bin
                        Integer experimentTreeNodeId = mapPPLNodes.get(jplacePhyloTreeNodeId);
                        //calculate the distance between these 2 nodeIds
                        //i.e. use the DTx and D'Tx matrices

                        //TODO: add the Nx ids in the expected placement binary
                        
                        //got coordinates of placed read
                        String[] readInfos=name.split("_");
                        int readStart=Integer.decode(readInfos[readInfos.length-2]);
                        int readEnd=Integer.decode(readInfos[readInfos.length-1]);

                        int nodeDistance = Dtx.getNodeDistance(prunedNodeId, experimentTreeNodeId);

                        //System.out.println(name+" -> nodeDistance:"+nodeDistance);
                        PPLResults.put(experimentLabel+":"+name, nodeDistance);
                        //software;Ax;k;alpha;Rx;read;node_dist
                        bw.append("PPL;"+experimentLabel+";;;;"+readLabel+";"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
                    }

                }
            }
            
            if (dg.doRAP) {

            
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
                    String kAlphaLabel=currentJPlaceFile.getParent().toFile().getName(); //1 times get parent  kx_ax/*.jplace
                    String[] data =kAlphaLabel.split("_");
                    int k=Integer.parseInt(data[0].substring(1));
                    float alpha=Float.parseFloat(data[1].substring(1));

                    //experiment
                    String experimentLabel=currentJPlaceFile.getParent().getParent().getFileName().toString(); //2 times get parent  Ax_nxx_xxx/kx_ax/logs/jplace
                    System.out.print("experimentLabel:"+experimentLabel);
                    int pruningNumber=Integer.parseInt(experimentLabel.split("_")[0].substring(1));
                    int prunedNodeId=Integer.parseInt(experimentLabel.split("_")[1].substring(2)); //i.e. Nx
                    //elements in filename in rappas ex: placements_ R3_nx3_la_r900.fasta_medium .jplace
                    String infos=currentJPlaceFile.getFileName().toString().split("\\.jplace$")[0].substring(11);  //11, to remove prefix "placements_"
                    String readLabel=infos.split("\\.")[0];
                    //System.out.print(" read:"+readLabel);
                    //elements in filename in rappas ex: R3 nx3 la r900.fasta medium
                    String[]elts=infos.split("_");
                    String dbSize=elts[elts.length-1];
                    //System.out.print(" dbSize:"+dbSize);
                    String readSizeElts=elts[elts.length-2];
                    String readSize=readSizeElts.substring(1, readSizeElts.length()-6);
                    //System.out.println(" readSize:"+readSize);


                    //pruning infos
                    int expectedPlacementIndex=NxIndex.get(prunedNodeId);

                    //load tree and expectedPlacement related to this Ax
                    PhyloTree experimentTree=experimentTrees.get(expectedPlacementIndex);
                    ArrayList<Integer> expectedPlacementNodeIds=expectedPlacementsNodeIds.get(expectedPlacementIndex);
                    //System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("experimentTree best placement(s):"+expectedPlacementNodeIds);
                    JplacerLoader RAPJplace=new JplacerLoader(currentJPlaceFile.toFile(), false);
                    //System.out.println("RAPJplace tree: "+RAPJplace.getTree());
                    //System.out.println("experimentTree: "+experimentTree);
                    if (RAPJplace.getTree().getNodeCount()!=experimentTree.getNodeCount()) {
                        System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                        System.out.println("They do not include the same trees for the same Nx experiment.");
                        System.exit(1);
                    }

                    //map jplace to experimentTree map(jplace nodeId)=experiment tree node id
                    HashMap<Integer, Integer> mapRAPNodes = RAPJplace.getTree().mapNodes(experimentTree);
                    //System.out.println("mapPPLNodes:"+mapPPLNodes);

                    //retrieve best placements
                    HashMap<String, ArrayList<Placement>> RAPBestPlacements = RAPJplace.getPlacements();
                    //System.out.println("RAPBestPlacements:"+RAPBestPlacements);

                    //for each placement (json item 'p' in the jplace)
                    for (Iterator<String> iterator = RAPBestPlacements.keySet().iterator(); iterator.hasNext();) {
                        String name = iterator.next();
                        //get best placement as the nodeId of the phylotree generated 
                        //during jplace parsing
                        Integer jplacePhyloTreeNodeId = RAPBestPlacements.get(name).get(0).getNodeId();
                        //verify if following placements are not same value
                        if (RAPBestPlacements.get(name).size()>1) {
                            if (RAPBestPlacements.get(name).get(1).getWeightRatio()==RAPBestPlacements.get(name).get(0).getWeightRatio()) {
                                System.out.println("!!!!!!!!!!!!!!  Identical weight ratios !");
                                System.out.println("Read: "+name);
                                System.out.println("File: "+currentJPlaceFile.toFile().getAbsolutePath());
                                System.exit(1);
                            }
                            
                        }
                        
                        //get its equivalent nodeId in the phylotree loaded from the 
                        //expected_placements.bin
                        Integer observedExperimentNodeId = mapRAPNodes.get(jplacePhyloTreeNodeId);
                        //calculate the distance between these 2 nodeIds
                        //i.e. use the DTx and D'Tx matrices
                        //TODO: add the Nx ids in the expected placement binary
                        
                        int nodeDistance = Dtx.getNodeDistance(prunedNodeId, observedExperimentNodeId);
                        
                        //test to determine if placement shifted to ancestors
                        //to do that, we check only distances of 1 from expected
                        //placement. 4 possible results: P,B,L,R
                        //
                        //        | P(arent)
                        //        |
                        //       / \\
                        //      /   \\  <--E(xpected placement)
                        //B(rother)  /\
                        //          /  \
                        //     L(eft)   R(ight)
                        //
                        PhyloNode observedPlacement = experimentTree.getById(observedExperimentNodeId);
                        PhyloNode expectedPlacement1 = experimentTree.getById(expectedPlacementNodeIds.get(0));
                        PhyloNode expectedPlacement2 =null;
                        if (expectedPlacementNodeIds.size()>1) {
                            expectedPlacement2 = experimentTree.getById(expectedPlacementNodeIds.get(1));
                        }
                        
                        boolean wasDistOne=false;
                        if (observedPlacement.getParent()!=null) {
                            //observed placement is expected placement
                            if ( expectedPlacementNodeIds.contains(observedExperimentNodeId) ) {
                                //case E
                                bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";E\n");
                                wasDistOne=true;
                            } else {
                                
                                //////////////
                                //check for expectedPlacement1
                                
                                //is not a leaf ? test if P
                                if(observedPlacement.getChildCount()>0) {
                                    Enumeration children = observedPlacement.children();
                                    while (children.hasMoreElements()) {
                                        PhyloNode nextElement = (PhyloNode)children.nextElement();
                                        //if a son of observed placement is expected placement, then observed is P                             
                                        if (nextElement.getId()==expectedPlacement1.getId()) {
                                            bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";P\n");
                                            wasDistOne=true;
                                            break;
                                        }
                                    }
                                } 
                                // if observed is son of expected, then L or R
                                if ( ((PhyloNode)observedPlacement.getParent()).getId()==expectedPlacement1.getId()) {
                                    if ( ((PhyloNode)observedPlacement.getParent().getChildAt(0)).getId()==observedPlacement.getId()) {
                                        bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";L\n");
                                        wasDistOne=true;
                                    } else {
                                        bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";R\n");
                                        wasDistOne=true;
                                    }
                                }
                                // test brothers
                                if ( ((PhyloNode)observedPlacement.getParent().getChildAt(0)).getId()==expectedPlacement1.getId()) {
                                    //observed is B on the right
                                    bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";BR\n");
                                    wasDistOne=true;
                                } else if ( ((PhyloNode)observedPlacement.getParent().getChildAt(1)).getId()==expectedPlacement1.getId() ) {
                                    bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";BL\n");
                                    wasDistOne=true;
                                }
                                
                                //////////////
                                //check for expectedPlacement2
                                
                                if (expectedPlacement2 !=null) {
                                    
                                    //is not a leaf ? test if P
                                    if(observedPlacement.getChildCount()>0) {
                                        Enumeration children = observedPlacement.children();
                                        while (children.hasMoreElements()) {
                                            PhyloNode nextElement = (PhyloNode)children.nextElement();
                                            //if a son of observed placement is expected placement, then observed is P                             
                                            if (nextElement.getId()==expectedPlacement2.getId()) {
                                                bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";P\n");
                                                wasDistOne=true;
                                                break;
                                            }
                                        }
                                    } 
                                    // if observed is son of expected, then L or R
                                    if ( ((PhyloNode)observedPlacement.getParent()).getId()==expectedPlacement2.getId()) {
                                        if ( ((PhyloNode)observedPlacement.getParent().getChildAt(0)).getId()==observedPlacement.getId()) {
                                            bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";L\n");
                                            wasDistOne=true;
                                        } else {
                                            bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";R\n");
                                            wasDistOne=true;
                                        }
                                    }
                                    // test brothers
                                    if ( ((PhyloNode)observedPlacement.getParent().getChildAt(0)).getId()==expectedPlacement2.getId() ) {
                                        //observed is B on the right
                                        bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";BR\n");
                                        wasDistOne=true;
                                    } else if ( ((PhyloNode)observedPlacement.getParent().getChildAt(1)).getId()==expectedPlacement2.getId() ) {
                                        bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";BL\n");
                                        wasDistOne=true;
                                    }
                                    
                                }
                            }
                            
                            
                            
                        } else {
                            //observed placement is root, case ignored
                        }
                        
                        if(!wasDistOne) {
                            bw3.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";N\n");
                        }
                        
                        //got coordinates of placed read
                        String[] readInfos=name.split("_");
                        int readStart=Integer.decode(readInfos[readInfos.length-2]);
                        int readEnd=Integer.decode(readInfos[readInfos.length-1]);
                        

                        //System.out.println(name+" -> nodeDistance:"+nodeDistance);
                        //software;Ax;k;alpha;Rx;read;node_dist
                        bw.append("RAP;"+experimentLabel+";"+k+";"+alpha+";"+dbSize+";"+readLabel+";"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");

                        
                        int distEPA=-1;
                        if (EPAResults.containsKey(experimentLabel+":"+name)) {
                            distEPA=EPAResults.get(experimentLabel+":"+name);
                        }
                        int distPPL=-1;
                        if (PPLResults.containsKey(experimentLabel+":"+name)) {
                            distPPL=PPLResults.get(experimentLabel+":"+name);
                        }
                        bw2.append(experimentLabel+";"+name+";"+k+";"+alpha+";"+dbSize+";"+readSize+";"+
                                    distEPA + ";" +
                                    distPPL + ";" +
                                    nodeDistance + "\n");
                    }

                } 
            
            }
            
            bw.close();
            bw2.close();
            bw3.close();
            
            System.out.println("DONE");
            
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
    
    
    public static void testPosteriorDFS(PhyloNode n) {
        Enumeration children = n.children();
        while (children.hasMoreElements()) {
            PhyloNode nextElement = (PhyloNode)children.nextElement();
            testPosteriorDFS(nextElement);
        }
        System.out.println(n);

    }
    
    public static void testAnteriorDFS(PhyloNode n) {
        System.out.println(n);
        Enumeration children = n.children();
        while (children.hasMoreElements()) {
            PhyloNode nextElement = (PhyloNode)children.nextElement();
            testPosteriorDFS(nextElement);
        }
    }
    
    
}
