
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

/**
 *
 * @author benjamin linard
 */
public class DistanceGenerator_LITE {
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME);

    //which component to score ?
    boolean doEPA=true;
    boolean doPPL=true;
    boolean doRAP=true;
    boolean doEPANG=true;
    boolean doAPPLES=true;
    
    //trifurcations follow naming convention :
    //Tx_trifuXX_nxXX_laXX
    boolean trifurcations=false;
    ArrayList<Integer> trifurcationsNxIndexes=null;
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose","1");
        
        HashMap<String,Integer> PPLResults=new HashMap<String,Integer>(); //map(experimentlabel_read_name)=ARTree_nodeid
        HashMap<String,Integer> EPAResults=new HashMap<String,Integer>(); //map(experimentlabel_read_name)=node_dist
        ObjectInputStream ois=null;
        
        try {
            
            System.out.println("ARGS: workDir doEPA[1/0] doEPANG[0/1] doPPL[1/0] doRAP[0/1] doAPPLES[0,1] [trifurcations:-1=no|1,45,48=list of Nx to test]");
            
            //launch
            DistanceGenerator_LITE dg=new DistanceGenerator_LITE();
            
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                dg.workDir=new File(args[0]);
                System.out.println("workDir: "+dg.workDir);
                if (Integer.parseInt(args[1])<1) { dg.doEPA=false; }
                if (Integer.parseInt(args[2])<1) { dg.doEPANG=false; }
                if (Integer.parseInt(args[3])<1) { dg.doPPL=false; }
                if (Integer.parseInt(args[4])<1) { dg.doRAP=false; }
                if (Integer.parseInt(args[5])<1) { dg.doAPPLES=false; }
                if (!args[6].equals("-1")) {
                    dg.trifurcations=true;
                    String[] trifuIndexes=args[6].split(",");
                    dg.trifurcationsNxIndexes=new ArrayList<>(trifuIndexes.length);
                    for (String trifuIndexe : trifuIndexes) {
                        dg.trifurcationsNxIndexes.add(Integer.valueOf(trifuIndexe));
                    }
                }
            }  
            System.out.println("trifurcations:"+dg.trifurcations);
            System.out.println("trifurcations Nx tested: "+dg.trifurcationsNxIndexes);
            
            //expected placement
            File expPLaceFile=new File(dg.workDir+File.separator+"expected_placements.bin");
            //Dtx
            File DtxFile=new File(dg.workDir+File.separator+"Dtx.csv");
           

            //load the expected placement
            ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(expPLaceFile),4096));
            System.out.println("Loading NxIndex");
            HashMap<Integer,Integer> NxIndex = (HashMap<Integer,Integer>)ois.readObject();
            System.out.println("Loading pruningIndex");
            HashMap<Integer,Integer> pruningIndex = (HashMap<Integer,Integer>)ois.readObject();
            System.out.println("Loading expected placements");
            ArrayList<ArrayList<Integer>> expectedPlacements = (ArrayList<ArrayList<Integer>>)ois.readObject();
            System.out.println("Loading trees");
            ArrayList<PhyloTree> experimentTrees = (ArrayList<PhyloTree>)ois.readObject();
            System.out.println("Loading trees trifurcations");
            ArrayList<ArrayList<PhyloTree>> experimentTreesTrifurcations = (ArrayList<ArrayList<PhyloTree>>)ois.readObject();      
            
            System.out.println("################################################");
            System.out.println("NxIndex="+NxIndex);
            System.out.println("pruningIndex="+NxIndex);
            System.out.println("expectedPlacements="+expectedPlacements);
            System.out.println("prunedTrees="+experimentTrees.size());
            for (int i = 0; i < experimentTrees.size(); i++) {
                PhyloTree get = experimentTrees.get(i);
                System.out.println(i+"th tree size test:"+get.getNodeCount());
            }
            if (dg.trifurcations) {
                for (int i = 0; i < dg.trifurcationsNxIndexes.size(); i++) {
                    int Nx=dg.trifurcationsNxIndexes.get(i);
                    System.out.println("Nx="+Nx+" with #trifurcationTrees="+experimentTreesTrifurcations.get(Nx).size());
                }
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
            if (dg.trifurcations) {
                bw.append("software;trifu;Ax;k;omega;read;readSize;node_dist;readStart;readEnd\n");
            } else {
                bw.append("software;Ax;k;omega;read;readSize;node_dist;readStart;readEnd\n");
            }
            
            
            
            BufferedWriter bwTrifu=null;
            BufferedWriter bw2Trifu=null;
            BufferedWriter bw3Trifu=null;
            if (dg.trifurcations) {
                //prepare a nice CSV file in which all data will be saved
                csvResult=Paths.get(dg.workDir.getAbsolutePath(),"results_trifu.csv");
                bwTrifu= Files.newBufferedWriter(csvResult);
                //header
                bwTrifu.append("software;Ax;k;omega;read;readSize;node_dist;readStart;readEnd\n");
               
            }
            
            
            
            
            if (dg.doEPA) {
                System.out.println("##############");
                System.out.println("## EPA");
                File EPAxDir=new File(dg.workDir+File.separator+"EPA");
                //load EPA jplace results
                List<Path> EPAJPlaceFiles = Files.find(EPAxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

                //for each jplace file, calculate node dist to expected placement
                for (int i = 0; i < EPAJPlaceFiles.size(); i++) {
                    Path currentJPlaceFile = EPAJPlaceFiles.get(i);
                    //SECOND element in filename in pplacer ex: RAxML_portableTree. R63_nx63_la1a.HCV_1_r200 .jplace
                    String jplaceLabel=currentJPlaceFile.getFileName().toString().split("\\.jplace$")[0]; 
                    String[]elts=jplaceLabel.split("_");
                    String readSize=elts[1].substring(1);
                    int pruning=Integer.valueOf(elts[0]);
                    System.out.println("--------------------------------------");
                    System.out.println("pruning:"+pruning+" read:"+readSize);

                    //load tree and expectedPlacements related to this Ax
                    PhyloTree experimentTree=experimentTrees.get(pruning);
                    experimentTree.initIndexes();
                    ArrayList<Integer> experimentPlacements=expectedPlacements.get(pruning);
                    //System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("experimentTree best placement(s):"+experimentPlacements);

                    
                    JplacerLoader EPAJplace=new JplacerLoader(currentJPlaceFile.toFile(), false);
                    
                    
                    if (EPAJplace.getTree().getNodeCount()!=experimentTree.getNodeCount()) {
                        EPAJplace.getTree().displayTree();
                        experimentTree.displayTree();
                        System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                        System.out.println("They do not include the same trees for the same Nx experiment.");
                        return;
                        //System.exit(1);
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
                        int experimentTreeNodeId = mapEPANodes.get(jplacePhyloTreeNodeId);
                        //System.out.println("experimentTreeNodeId: "+experimentTreeNodeId);
                        //calculate the distance between these 2 nodeIds
                        //i.e. use the DTx and D'Tx matrices
                        //System.out.println("expected:"+expectedPlacements.get(pruning));
                        
                        int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), experimentTreeNodeId);
                        
                        //got coordinates of placed read
                        String[] readInfos=name.split("_");
                        long readStart=0;
                        long readEnd=0;
                        try {
                            readStart=Long.decode(readInfos[readInfos.length-2]);
                        } catch (NumberFormatException ex) {
                            ex.printStackTrace();
                        }
                        try  {
                            readEnd=Long.decode(readInfos[readInfos.length-1]);
                        } catch (NumberFormatException ex) {
                            ex.printStackTrace();
                        }

                        //System.out.println(name+" -> nodeDistance:"+nodeDistance);
                        EPAResults.put(pruning+":"+name, nodeDistance);
                        //software;Ax;k;omega;Rx;read;node_dist
                        bw.append("EPA;"+pruning+";;;"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
                    }

                }
                bw.flush();
            }    
            
            if (dg.doEPANG) {
                System.out.println("##############");
                System.out.println("## EPA-ng");
                File EPANGxDir=new File(dg.workDir+File.separator+"EPANG");
                //load EPA jplace results
                List<Path> EPANGJPlaceFiles = Files.find(EPANGxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

                //for each jplace file, calculate node dist to expected placement
                for (int i = 0; i < EPANGJPlaceFiles.size(); i++) {
                    Path currentJPlaceFile = EPANGJPlaceFiles.get(i);
                    //FIRST  element in filename in epa-ng ex: R0_nx3_la95_r150. aln. fasta_queries_epa_result. jplace
                    String jplaceLabel=currentJPlaceFile.getFileName().toString().split("\\.jplace$")[0]; 
                    String[]elts=jplaceLabel.split("_");
                    String readSize=elts[1].substring(1);
                    int pruning=Integer.valueOf(elts[0]);
                    System.out.println("--------------------------------------");
                    System.out.println("pruning:"+pruning+" read:"+readSize);

                    //load tree and expectedPlacements related to this Ax
                    PhyloTree experimentTree=experimentTrees.get(pruning);
                    experimentTree.initIndexes();
                    ArrayList<Integer> experimentPlacements=expectedPlacements.get(pruning);
                    //System.out.println("experimentTree nodeIds:"+experimentTree.getNodeIdsByDFS());
                    //System.out.println("experimentTree best placement(s):"+expectedPlacementNodeIds);

                    
                    JplacerLoader EPANGJplace=new JplacerLoader(currentJPlaceFile.toFile(), false);

                    //version of EPA-ng prior to 0.3.4 unroots the input tree
                    //this needs to be corrected, posterior versions, not need to correct                    
                    
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
                        int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), experimentTreeNodeId);
                        
                        //got coordinates of placed read
                        String[] readInfos=name.split("_");
                        long readStart=0;
                        long readEnd=0;
                        try {
                            readStart=Long.decode(readInfos[readInfos.length-2]);
                        } catch (NumberFormatException ex) {
                            ex.printStackTrace();
                        }
                        try  {
                            readEnd=Long.decode(readInfos[readInfos.length-1]);
                        } catch (NumberFormatException ex) {
                            ex.printStackTrace();
                        }

                        //System.out.println(name+" -> nodeDistance:"+nodeDistance);
                        EPAResults.put(pruning+":"+name, nodeDistance);
                        //software;Ax;k;omega;Rx;read;node_dist
                        bw.append("EPA-ng;"+pruning+";;;"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
                    }

                }
                bw.flush();
            }    
            
            
            if (dg.doPPL) {
                
                System.out.println("##############");
                System.out.println("## PPL");
                File PPLxDir=new File(dg.workDir+File.separator+"PPLACER");
                //load EPA jplace results
                List<Path> PPLJPlaceFiles = Files.find(PPLxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

                //for each jplace file, calculate node dist to expected placement
                for (int i = 0; i < PPLJPlaceFiles.size(); i++) {
                    Path currentJPlaceFile = PPLJPlaceFiles.get(i);
                    String jplaceLabel=currentJPlaceFile.getFileName().toString().split("\\.jplace$")[0]; 
                    String[]elts=jplaceLabel.split("_");
                    String readSize=elts[1].substring(1);
                    int pruning=Integer.valueOf(elts[0]);
                    System.out.println("--------------------------------------");
                    System.out.println("pruning:"+pruning+" read:"+readSize);

                    //load tree and expectedPlacements related to this Ax
                    PhyloTree experimentTree=experimentTrees.get(pruning);
                    ArrayList<Integer> expectedPlacementNodeIds=expectedPlacements.get(pruning);
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
                        long readStart=0;
                        long readEnd=0;
                        try {
                            readStart=Long.decode(readInfos[readInfos.length-2]);
                        } catch (NumberFormatException ex) {
                            ex.printStackTrace();
                        }
                        try  {
                            readEnd=Long.decode(readInfos[readInfos.length-1]);
                        } catch (NumberFormatException ex) {
                            ex.printStackTrace();
                        }

                        int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), experimentTreeNodeId);

                        //System.out.println(name+" -> nodeDistance:"+nodeDistance);
                        PPLResults.put(pruning+":"+name, nodeDistance);
                        //software;Ax;k;omega;Rx;read;node_dist
                        bw.append("PPL;"+pruning+";;;"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
                    }

                }
                bw.flush();
            }
            
            if (dg.doRAP) {

            
                System.out.println("##############");
                System.out.println("## RAP");
                File DxDir=new File(dg.workDir+File.separator+"RAPPAS");
                //load EPA jplace results
                ArrayList<File> RAPResults=new ArrayList<>();
                List<Path> RAPPJPlaceFiles = Files.find(DxDir.toPath(), 999, (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")).collect(Collectors.toList());

                //for each jplace file, calculate node dist to expected placement
                for (int i = 0; i < RAPPJPlaceFiles.size(); i++) {
                    Path currentJPlaceFile = RAPPJPlaceFiles.get(i);
                    System.out.println(currentJPlaceFile.toAbsolutePath());
                    //not related to trifurcations
                    if (!currentJPlaceFile.getFileName().toString().contains("trifu_")) {
                        buildRAPPASNodeDistances(false, currentJPlaceFile,Dtx,NxIndex,expectedPlacements,experimentTrees,experimentTreesTrifurcations,EPAResults,PPLResults,bw,pruningIndex);
                    } else {//related to trifurcations
                        if (dg.trifurcations) {
                            buildRAPPASNodeDistances(true, currentJPlaceFile,Dtx,NxIndex,expectedPlacements,experimentTrees,experimentTreesTrifurcations,EPAResults,PPLResults,bwTrifu,pruningIndex);
                        }
                    }
                } 
            
            }
            
            bw.close();
            
            System.out.println("DONE");
            
            System.exit(0);
            
            
        
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DistanceGenerator_LITE.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException | ClassNotFoundException ex) {
            Logger.getLogger(DistanceGenerator_LITE.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                ois.close();
            } catch (IOException ex) {
                Logger.getLogger(DistanceGenerator_LITE.class.getName()).log(Level.SEVERE, null, ex);
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
    
    /**
     * operations related to RAPPAS jplace parsing, node distance, 
     * and comparisons to EPA / PPlacer
     * @param currentJPlaceFile
     * @param Dtx
     * @param NxIndex
     * @param expectedPlacementsNodeIds
     * @param experimentTrees
     * @param experimentTreesTrifurcations
     * @param EPAResults
     * @param PPLResults
     * @param bw
     * @param bw2
     * @param bw3
     * @throws IOException 
     */
    private static void buildRAPPASNodeDistances(  
                                            boolean trifurcations,
                                            Path currentJPlaceFile, 
                                            Dtx Dtx, 
                                            HashMap<Integer,Integer> NxIndex,
                                            ArrayList<ArrayList<Integer>> expectedPlacementsNodeIds,
                                            ArrayList<PhyloTree> experimentTrees,
                                            ArrayList<ArrayList<PhyloTree>> experimentTreesTrifurcations,
                                            HashMap<String,Integer> EPAResults,
                                            HashMap<String,Integer> PPLResults,
                                            BufferedWriter bw, 
                                            HashMap<Integer,Integer> pruningIndex
            
                                            ) throws IOException {
        
        String jplaceLabel=currentJPlaceFile.getFileName().toString(); //2 times get parent  Ax_nxx_xxx/kx_ax/logs/jplace
        System.out.print("jplaceLabel:"+jplaceLabel);
        String[] data =currentJPlaceFile.toFile().getName().split("_");
        
        //k and omega
        int k=Integer.parseInt(data[2].substring(1));
        float omega=Float.parseFloat(data[3].substring(1));
        //experiment
        int pruning=Integer.parseInt(data[0]);
        String readSize=data[1].substring(1);
        //elements in filename in rappas ex: placements_ R3_nx3_la_r900.fasta_medium .jplace
        String infos=currentJPlaceFile.getFileName().toString().split("\\.jplace$")[0].substring(11);  //11, to remove prefix "placements_"
        String readLabel=infos.split("\\.")[0];
        System.out.println("pruning: "+pruning);
        System.out.println("readSize: "+readSize);
        System.out.println("k: "+k);
        System.out.println("omega: "+omega);


        //pruning infos

        //load tree and expectedPlacements related to this Ax
        PhyloTree experimentTree=experimentTrees.get(pruning);
        int trifuNumber=-1;
        if (trifurcations) {
            //trifu_X_
            trifuNumber=Integer.parseInt(currentJPlaceFile.getFileName().toString().split("_")[1]);
            experimentTree=experimentTreesTrifurcations.get(pruning).get(trifuNumber);
            System.out.println("ND for trifu: "+trifuNumber+"  root: "+experimentTree.getRoot().toString());
        }
        
        
        ArrayList<Integer> expectedPlacementNodeIds=expectedPlacementsNodeIds.get(pruning);
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
                    //System.exit(1);
                }
            }

            //get its equivalent nodeId in the phylotree loaded from the 
            //expected_placements.bin
            Integer observedExperimentNodeId = mapRAPNodes.get(jplacePhyloTreeNodeId);
            //calculate the distance between these 2 nodeIds
            //i.e. use the DTx and D'Tx matrices
            //TODO: add the Nx ids in the expected placement binary

            int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), observedExperimentNodeId);

        

            //got coordinates of placed read
            String[] readInfos=name.split("_");
            long readStart=0;
            long readEnd=0;
            try {
                readStart=Long.decode(readInfos[readInfos.length-2]);
            } catch (NumberFormatException ex) {
                ex.printStackTrace();
            }
            try  {
                readEnd=Long.decode(readInfos[readInfos.length-1]);
            } catch (NumberFormatException ex) {
                ex.printStackTrace();
            }


            //System.out.println(name+" -> nodeDistance:"+nodeDistance);
            //software;Ax;k;omega;Rx;read;node_dist
            if (!trifurcations) {
                bw.append("RAP;"+pruning+";"+k+";"+omega+";"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
            } else {
                bw.append("RAP;"+trifuNumber+";"+pruning+";"+k+";"+omega+";"+name+";"+readSize+";"+nodeDistance+";"+readStart+";"+readEnd+"\n");
            }


        }

        
    }
    
    
}
