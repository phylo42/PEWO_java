
import alignement.Alignment;
import core.AAStates;
import core.DNAStatesShifted;
import core.States;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.file.attribute.PosixFilePermission;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.PhyloTree.Path;
import tree.PhyloTreeModel;

import javax.swing.tree.TreeNode;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * From an alignment and tree, build directories of pruning experiments, 
 if tree is unrooted then roots it for each pruning experiment,
 saves pruned alignments in AxFile, pruned trees in TxFile and generated reads in Rx,
 build Dtx matrice which defines node_dist between edges (note that roots 
 are not considered in this distance, the 2 edges sons of root are considered
 as the same placement) 
 * @author ben
 */
public class PrunedTreeGenerator_LITE {
    
    //pseudorandom generator
    Long seed=1L;
    Random rand = new Random(seed);
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME);
    
    //default dataset
    File alignFile=new File(HOME);
    File treeFile=new File(HOME);
    
    //States
    States s=new DNAStatesShifted();
    
    //nodes selected for pruning at this launch
    Integer[] prunedNodeIds=null;
    //nodes effectively pruned because didn'TxFile produce tre of less than 4 leaves
    ArrayList<Integer> actuallyPrunedNodeIdsIndexes=new ArrayList<>();
    
    //list of new files
    public HashMap<Integer,List<String>> readFiles = new HashMap<>(); //map of reads, map(0)=[A0_nx4_la_r150,A0_nx4_la_r300]
    public ArrayList<PhyloTree> prunedTrees = new ArrayList<>(); //list of pruned Trees
    public ArrayList<ArrayList<PhyloTree>> prunedTreesTrifurcations = new ArrayList<>(); //list of pruned Trees
    public File fileDtx=null;
    public File fileD2tx=null;
    
    //pruning fraction
    //double percentPruning=0.5; //10%
    int pruningCount=100;

    //read generation: nomral distrib around mean R with sd (R/4)
    //and min length Rmin
    int[] R={3,6};
    //standard deviation for each r in R
    double Rsd=0.5;
    int Rmin=75; //let's consider that we have at least 75bp reads
    
    //set if analysis is protein or DNA/RNA
    boolean proteinAnalysis=false;
    
    //set the critera for branch injection
    float minBranchLength=-1.0f; //basically all branches, of all length (even 0)
    int branchPerEdge=1;


    
    //file permissions given to the qsub scripts
    EnumSet<PosixFilePermission> perms =
            EnumSet.of(
                        PosixFilePermission.OWNER_READ,
                        PosixFilePermission.OWNER_WRITE,
                        PosixFilePermission.OWNER_EXECUTE,
                        PosixFilePermission.GROUP_READ,
                        PosixFilePermission.GROUP_EXECUTE,
                        PosixFilePermission.OTHERS_READ,
                        PosixFilePermission.OTHERS_EXECUTE
                    );
        
    public static void main(String[] args) {
        
        System.out.println("ARGS: workDir align tree #prunings(int) readSize1(int),readSize2(int),... readSD(int) branchPerEdge(int) [nucl=0|prot=1] ");
        
        System.out.println("Command: "+Arrays.toString(args).replaceAll(",", " "));

        
        try {
            //launch
            PrunedTreeGenerator_LITE ptg=new PrunedTreeGenerator_LITE();
            
            //PREPARE DATA SOURCES
            ///////////////////////////////////////////////////
            if (args.length>0) {
                ptg.workDir=new File(args[0]);
                ptg.alignFile=new File(args[1]);
                ptg.treeFile=new File(args[2]);
                ptg.pruningCount=Integer.parseInt(args[3]);
                String[] readSizes=args[4].split(",");
                ptg.R=new int[readSizes.length];
                for (int i = 0; i < readSizes.length; i++) {
                    ptg.R[i]=Integer.valueOf(readSizes[i]);
                }
                ptg.Rsd=Integer.parseInt(args[5]);
                //optionnals, if no args use default values
                ptg.branchPerEdge=Integer.parseInt(args[6]);
                int protein=Integer.parseInt(args[7]);
                ptg.proteinAnalysis=(protein>0);

                
            }
            
            
            
            System.out.println("workDir: "+ptg.workDir);
            System.out.println("alignFile: "+ptg.alignFile);
            System.out.println("treeFile: "+ptg.treeFile);
            System.out.println("pruningCount: "+ptg.pruningCount);
            System.out.println("readSizes: "+Arrays.toString(ptg.R));
            System.out.println("branchPerEdge: "+ptg.branchPerEdge);
            System.out.println("protein:"+ptg.proteinAnalysis);

     
            //LOAD TREE / ALIGNMENTS
            ///////////////////////////////////////////////////
            //load alignment
            ptg.s=new DNAStatesShifted();
            if (ptg.proteinAnalysis) {
                ptg.s=new AAStates(true);
            }
            
            FASTAPointer fp=new FASTAPointer(ptg.alignFile, false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
            }
            Alignment align=new Alignment(ptg.s,fastas);
            System.out.println(align.describeAlignment(false));
            fp.closePointer();
            //load ancTree
            BufferedReader br=new BufferedReader(new FileReader(ptg.treeFile));
            String line=null;
            String treeString=null;
            while((line=br.readLine())!=null) {treeString=line;}
            System.out.println("Force tree rooting (if unrooted)");
            PhyloTree tree = NewickReader.parseNewickTree2(treeString, true, false);

            tree.initIndexes();
            System.out.println("Original tree, # nodes: "+tree.getNodeCount());
            System.out.println("Is rooted: "+tree.isRooted());
            //tree.displayTree();

            //test if alignment/tree labels are matching.
            //if not, exit before raising more errors later in the algo...
            List<String> alignLabels = Arrays.asList(align.getRowLabels());
            int notFoundCount=0;
            for (Iterator<String> iterator = alignLabels.iterator(); iterator.hasNext();) {
                String next = iterator.next();
                if (!tree.getLabelsByDFS().contains(next)) {
                    System.out.println("Alignment label \""+next+"\" not found in given tree labels.");
                    notFoundCount++;
                }
            }
            if (notFoundCount>0) {System.exit(1);}
            alignLabels=null;
            
            //also tests if some nodes got erroneously parsed because the newick
            //tree was not valid. This results to "null" labels
            ArrayList<Integer> nodeIdsByDFS = tree.getNodeIdsByDFS();
            for (int id:nodeIdsByDFS) {
                if (tree.getById(id).getLabel()==null) {
                    System.out.println("Error encountered on tree node:"+tree.getById(id).toString());
                    System.out.println("Please check the integrity of your newick tree.");
                    System.exit(1);
                }
            }
            
            
            
            //PREPARE ALL PRUNING EXPERIMENT FILES (AxFile, TxFile, Rx)
            //A: pruned alignments
            //T: pruned trees
            //G: pruned leaves, as complete sequences
            //R: reads built from the pruned leaves
            ///////////////////////////////////////////////////
            if (!ptg.workDir.exists()) {
                ptg.workDir.mkdir();
            }
            ptg.generatePrunedTrees(ptg.workDir,tree,align);
            System.out.println("PRUNED TREE GENERATION DONE !");
            
            ////////////////////////////////////////////////////////////////////
            //build HMM alignments for this current pruning experiment
            //i.e. : build the hmm profile from the corresponding AxFile
            //       align all Rx reads to this profile
            //       convert the result from stockolm to fasta alignment

            
            
            
            Thread.sleep(1);
            
        } catch (IOException ex) {
            Logger.getLogger(PrunedTreeGenerator_LITE.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(PrunedTreeGenerator_LITE.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    
    /**
     * main script
     */
    private void generatePrunedTrees(File workDir,PhyloTree tree,Alignment align) throws IOException {
        
        //some definitions:
        //several nodeIds are selected randomly:                                prunedNodeIds=               [0,1,2,3,4,5]
        //each node Nx of this list will attempt a pruning experiment:          nx_nodeId=...
        //here, nx_nodeId=0 targets the root, which cannot be pruned,
        //this Nx pruning is SKIPPED
        //list of prunedNodeIds indexes that were actually pruned               actuallyPrunedNodeIdsIndexes=[1, 2, 3, 4, 5]
        //map(nodeId)=index                                                     NxIndex=                     {1=0, 2=1, 3=2, 4=3, 5=4}
        //expected placments (2 edges if non-rooted input was --froot)          expectedPlacements=          [[19, 48], [13], [10], [5], [4]]
        //PhyloTree list, saved to have same ids as in expectedPlacements       prunedTrees=                 [PhyloTree1,PhyloTree2,...,PhyloTree5]
        
        //Are saved in expected_placements.bin :
        // - NxIndex
        // - expectedPlacements
        // - prunedTrees        
        
        //==> modif pr que prunedAlignmentsFiles et prunedTreesFiles aient même taille que actuallyPrunedNodeIdsIndexes
        
        
        //states
        States states=new DNAStatesShifted();
        if (proteinAnalysis) {
            states=new AAStates(true);
        }
        
        //to ensure consistent direcory names (omega-> _ax.xx_)
        NumberFormat nf = NumberFormat.getNumberInstance(Locale.UK);
        nf.setMinimumFractionDigits(2);
        nf.setMaximumFractionDigits(2);
        
        //get list of internal Nx ids
        Integer[] nodeIds=new Integer[tree.getNodeIdsByDFS().size()];
        System.out.println("This tree contains "+nodeIds.length+" nodeIds.");
        //shuffle their order
        shuffleArray(nodeIds);
        //not shuffled
        //nodeIds=ancTree.getNodeIdsByDFS().toArray(nodeIds);
        //define first x% as pruning experiments
        if(pruningCount>nodeIds.length) {
            pruningCount=nodeIds.length;
        }
        prunedNodeIds=Arrays.copyOfRange(nodeIds, 0, pruningCount);
        System.out.println("# pruning attempts: "+prunedNodeIds.length);
        System.out.println("prunedNodeIds, selected for pruning attempts (raw): "+Arrays.toString(prunedNodeIds));
        //sort to make more comprehensive output matrices Dtx and D'tx
        Arrays.sort(prunedNodeIds);
        System.out.println("prunedNodeIds, selected for pruning attempts (ordered): "+Arrays.toString(prunedNodeIds));

        //log registering all skipped Nx, and the reason
        File skippedLog=new File(workDir+File.separator+"SKIPPED_Nx");
        BufferedWriter skippedLogBw=new BufferedWriter(new FileWriter(skippedLog,false));
        
        //prepare the output directories
        File AxDir=new File(workDir+File.separator+"A");
        File TxDir=new File(workDir+File.separator+"T");
        File RxDir=new File(workDir+File.separator+"R");
        File GxDir=new File(workDir+File.separator+"G");
        AxDir.mkdir();
        TxDir.mkdir();
        RxDir.mkdir();
        GxDir.mkdir();
            
        //write in a binary file the pruned ancTree and the expected placement,
        //that is the branch b_new (and b_new_p if rerooting), see algo below.
        //expected placement is saved through an integer array of 1 or 2 elements
        //integer is the n of the node son of b_new
        File expectedPlacementsFile=new File(workDir.getAbsolutePath()+File.separator+"expected_placements.bin");
        ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(expectedPlacementsFile),4096));
        HashMap<Integer,Integer> NxIndex=new HashMap<>(); //map(120,it.e. nx120)=0; map(1142)=1 ; map(5454)=2 ...
        HashMap<Integer,Integer> pruningIndex=new HashMap<>(); //map(0)=120,it.e. nx120; map(1)=1142 ; map(2)=5454 ...
        ArrayList<ArrayList<Integer>> expectedPlacements=new ArrayList<>(); //tab[1]=[124,123] ; tab[2]=[1143] ...
        
        //prepare Dtx and D'tx matrices, build their headers
        fileDtx=new File(workDir+File.separator+"Dtx.csv");
        fileD2tx=new File(workDir+File.separator+"D2tx.csv");
        BufferedWriter brDtx=new BufferedWriter(new FileWriter(fileDtx));
        BufferedWriter brD2tx=new BufferedWriter(new FileWriter(fileD2tx));
        StringBuilder nodeDistMatrix=new StringBuilder();
        StringBuilder branchDistMatrix=new StringBuilder();
        nodeDistMatrix.append("nodeLabels;");
        for (int n=0;n<nodeIds.length;n++)
            nodeDistMatrix.append(";"+tree.getById(nodeIds[n]).getLabel());
        nodeDistMatrix.append("\n");
        nodeDistMatrix.append(";nodeIds");
        for (int n=0;n<nodeIds.length;n++)
            nodeDistMatrix.append(";"+nodeIds[n]);
        nodeDistMatrix.append("\n");
        branchDistMatrix.append(new String(nodeDistMatrix)); //simple contructor copy
        

            
        //launch prunings for each selected Nx
        //////////////////////////////////////////////////////////////
        
        for (int i = 0; i < prunedNodeIds.length; i++) {
            Integer nx_nodeId = prunedNodeIds[i];
            //System.out.println("--------------------------------------");
            //System.out.println("copying ancTree and alignment before pruning...");
            PhyloNode rootCopy=tree.getRoot().copy();
            PhyloTree treeCopy=new PhyloTree(new PhyloTreeModel(rootCopy),tree.isRooted(), false);
            //System.out.println("indexing ancTree ...");
            treeCopy.initIndexes();
            //copying alignment
            Alignment alignCopy=align.copy();
            //System.out.println("Starting pruning...");
            
            //current root Nx defining the pruned clade
            PhyloNode Nx= treeCopy.getById(nx_nodeId);
            System.out.println("--------------------------------------");
            System.out.println("selected Nx: x="+i+"  node:"+Nx);
            if (Nx.isRoot()) {
                skippedLogBw.append("Nx="+i+"\t"+Nx.toString()+"\tIs root, so skipped.\n");
                System.out.println("SKIPPED: this node is root (id="+nx_nodeId+"), no pruning.");
                continue;
            }
                        
            
            ///////////////////////////////////////////////////////////
            //first, modify alignment and deleted subtree related to this Nx
            
            
            //if leaf, removes from multiple align
            ArrayList<String> leavesRemoved=new ArrayList<>();
            if (Nx.isLeaf()) {
                leavesRemoved.add(Nx.getLabel());
                alignCopy.removeSequence(Nx.getLabel());
            //if an internal Nx
            } else {
                //enumerate nodes in subtree
                PhyloNode nextNx =null;
                Enumeration<TreeNode> DFSenum = Nx.depthFirstEnumeration();
                //careful, ancTree topology change are reflected in the enumeration
                //so we need to remove node AFTER this this transversal postorder
                while (DFSenum.hasMoreElements()) {
                    nextNx=(PhyloNode)DFSenum.nextElement();
                    //if leaf, removes
                    if (nextNx.isLeaf()) {
                        leavesRemoved.add(nextNx.getLabel());
                        alignCopy.removeSequence(nextNx.getLabel());
                    } 
                }
                //delete Nx children (this does not delete the sub nodes from memory !!!)
                nextNx.removeAllChildren();
            }
            //if this removal ends to AxFile multiple alignment with only <3 leaves
            //we don'TxFile go further and pass to the next Nx.
            //This can happen for instance in this nx case:
            //
            //   R_____Leaf             or       R______Leaf1        
            //    |                               |   |_Leaf2
            //    |____Nx(subtree)                |_____Nx(subtree)
            //
            //
            //
            if (alignCopy.getRowLabels().length<3) {
                skippedLogBw.append("Nx="+i+"\t"+Nx.toString()+"\tPruning resulting to less than 3 leaves, so skipped.\n");
                System.out.println("SKIPPED: This pruning results to a 3 with less than 3 leaves !  --> Nx skipped");              
                continue;
            }
            System.out.println("LeavesRemoved: "+leavesRemoved.size());
            
            
            
            //REGISTER THIS PRUNING
            //if reach here, we made AxFile ancTree with more than 3 leaves, this pruning
            // will be effectively operated
            actuallyPrunedNodeIdsIndexes.add(i);
            //so let's register this pruning index in the expected_placement.bin file
            NxIndex.put(nx_nodeId,actuallyPrunedNodeIdsIndexes.size()-1);
            pruningIndex.put(actuallyPrunedNodeIdsIndexes.size()-1, nx_nodeId);
            

            
            File AxFile=new File(AxDir+File.separator+i+".align");
            File TxFile=new File(TxDir+File.separator+i+".tree");
            File GxFile=new File(GxDir+File.separator+i+".fasta");
            
            
            
            ////////////////////////////////////////////////////////////////////
            //second, built TxFile.
            //removal of node Nx of TxFile is operated on the basis of the following
            //definitions:
            // -Np    parent of Nx
            // -Np''  parent of Np
            // -Np'   son of Np
            //The goal: best placement of all Nx subclade leaves is the edge 
            //Np'-Np'', accessed through PhyloNode Np'.
            //      _          
            //     \ /                
            //      Np''               
            //      |                   
            //      Np         
            //     /  \                 
            //    Nx  Np'                
            //   /_\  /_\                  
            //
            //remove selected Nx itself and link its sister to parent
            PhyloNode Np=(PhyloNode)Nx.getParent(); //Np
            //System.out.println("Np:"+Np);
            float Np_bl=Np.getBranchLengthToAncestor(); //bl of Np to Np''
            PhyloNode Np_pp=(PhyloNode)Np.getParent(); //Np''
            //System.out.println("Np_pp "+Np_pp);
            PhyloNode Np_p=(PhyloNode)Np.getChildBefore(Nx); //Np' if on left
            if (Np_p==null) //if not on left, will be on right
                Np_p=(PhyloNode)Np.getChildAfter(Nx); //Np' if on right
            //System.out.println("Np_p "+Np_p);
            float Np_p_bl=Np_p.getBranchLengthToAncestor(); //bl of Np' to Np

            //memorize which edge is the best placement, this is AxFile node, 
            //but the edge going to its parent is best placement
            //if root is best placment, this edge is in fact null and of 
            //branch length =0
            PhyloNode b_new=null;
            PhyloNode b_new_p=null;
            
            //when Nx is son of root, we take the second son as Np
            //basically, in this case we just remove the subtree of Nx and
            // make Np the root; 
            //the best placement should then be the Np'-Np'' edge, which is 
            //split by the root, so both edge are best placement
            //                                 
            //                       
            //  added_root   ==>  Nx----Np       ==>      added_root
            //     /  \     shift      /  \                 /   \
            //    Nx   Np            Np'   Np''           Np'    Np''
            //   /_\  /_\           /_\    /_\           /_\     /_\
            //
            //
            if (Np_pp==null) { 
               
                //shift the Np pointer to the second root son
                if (Nx==Np.getChildAt(0)) {
                    Np=Np.getChildAt(1);
                } else {
                    Np=Np.getChildAt(0);
                }
                //set which is new Np' and Np''
                Np_p=Np.getChildAt(0);
                Np_pp=Np.getChildAt(1);
                //disconnect Nx and Np from root
                Nx.removeFromParent();
                Np.removeFromParent();
                //change Np as new root
                Np.setLabel("added_root");
                Np.setBranchLengthToAncestor(0.0f);
                //rebuilt phylotree using this new root
                treeCopy=new PhyloTree(new PhyloTreeModel(Np),true, false);
                //memorize best placement
                b_new=Np_p;
                b_new_p=Np_pp;
                

                
            //in all other case, the ancTree is disconnected and reconnected
            //      _                   _
            //     \ /                 \ /
            //      Np''                Np''
            //      |                   |
            //      Np          ==>     |
            //     /  \                  \
            //    Nx  Np'                 Np'
            //   /_\  /_\                 /_\    
            //
            } else{ 
                //disconnect all
                Np.removeFromParent();
                Np_p.removeFromParent();            
                Nx.removeFromParent();
                //connect Np' tp Np'' and update bl, this link is b_new
                Np_pp.add(Np_p);
                Np_p.setBranchLengthToAncestor(Np_bl+Np_p_bl);
                b_new=Np_p;
                b_new_p=null;
                //detach and get rid of Nx and Np
                Np.removeFromParent();
                Np=null;
            }
            
            //SAVE THE MODIFIED TREES AND ALIGNMENTS
            
            //write alignment to file and list
            //remove columns of 100% gaps only to not interfere in
            //ancestral reconstruction
            alignCopy.reduceAlignment(1.0);
            alignCopy.writeAlignmentAsFasta(AxFile);
            System.out.println("Write: "+AxFile.getAbsolutePath());
            //save the pruned ancTree
            System.out.println("Indexing pruned tree");
            treeCopy.initIndexes(); //necessary to remove former nodes from the maps shortcuts
            System.out.println("pruned tree(treeCopy), #nodes :"+treeCopy.getNodeCount());
            System.out.println("pruned tree(treeCopy), #leaves:"+treeCopy.getLeavesCount());
            //it is necessary to keep the trees objects, if saved as
            //newick, then reloaded, the ancTree node ids will be different
            //and the expectedPlacement map will not match
            prunedTrees.add(treeCopy);
            System.out.println("Write: "+TxFile.getAbsolutePath());
            NewickWriter nw=new NewickWriter(TxFile);            
            nw.writeNewickTree(treeCopy, true, true, false, false);  //no internal node names if PAML !
            nw.close();
            
            //SAVE THE EXPECTED PLACEMENT 
            //which is b_new (and b_new_p if rerooting)
            //fill the integer list
            ArrayList<Integer> array=new ArrayList<>();
            array.add(b_new.getId());
            if (b_new_p!=null)
                array.add(b_new_p.getId());
            expectedPlacements.add(array);
            
            
            ////////////////////////////////////////////////////////////////////
            //third, build Dtx and D'tx line corresponding to this pruning
            //note: for very big trees, should build that as an object
            //and save it by serialization ?
        
            nodeDistMatrix.append(Nx.getLabel()+";"+Nx.getId());
            branchDistMatrix.append(Nx.getLabel()+";"+Nx.getId());  //TODO: this still contains errors (dist X0 to neighboor node not taken into account for node)... use nodeDist for now
            //System.out.println(Arrays.toString(prunedNodeIds));
            //System.out.println("Np_p:"+Np_p);
            //System.out.println("Np_pp:"+Np_pp);
            //System.out.println("b_new:"+b_new);
            //using the nodeIds table, Dtx column order will match the
            //shuffled nodeIds, like this the first xx% correspond to the
            //pruned nodes. we could have used ancTree.getNodeIdsByDFS() too.
            for (int n=0;n<nodeIds.length;n++) {
                PhyloNode currentNode=treeCopy.getById(nodeIds[n]);
                //System.out.println("+++++++++++currentNode:"+currentNode);
                nodeDistMatrix.append(";");
                branchDistMatrix.append(";");
                if (currentNode==null) { //this node was pruned, so absent from TxFile tree copy
                    nodeDistMatrix.append("-1");
                    branchDistMatrix.append("-1.0");
                    continue;
                }  else if (currentNode==b_new) {//this is b_new itself, all dist=0
                    //System.out.println("b_new itself, dist=0");
                    nodeDistMatrix.append("0");
                    branchDistMatrix.append("0.0");
                    continue;
                } else {
                    //retrieve path from this new edge to all other nodes
                    Path shortestPath = treeCopy.shortestPath(treeCopy.getRoot(), b_new,currentNode);
                    //System.out.println("path: "+shortestPath);
                    
                    //corrections to node distance
                    int correctedNodeDistance=shortestPath.nodeDistance;
                    
                    //1st correction to node distance:
                    //added_root comes from forced rooting of the input ancTree
                    //and is injected betwen Np' and Np''
                    //if this triplet is on the path, a correction of -1
                    //has to be brought to nodeDistance (like if this added_root
                    //didn't existed)
                    if ( shortestPath.isWithAddedRoot() && !currentNode.getLabel().equals("added_root")) {
                        //System.out.println("ADDED_ROOT CORRECTION !");
                        correctedNodeDistance-=1;
                    }
                    //2nd correction to node distance:
                    //probably useless as added_root should be ignored in 
                    //distance matrices, but just in case
                    if ( currentNode.getLabel().equals("added_root")  ) {
                        //System.out.println("ADDED_ROOT ITSELF CORRECTION !");
                        correctedNodeDistance-=1;
                    }
                    //System.out.println("corrected: "+correctedNodeDistance);
                    nodeDistMatrix.append(correctedNodeDistance);
                    
                    
                    //for the branch ditance, we need to add the distance that
                    //was separating Np and Np'' if the path start to go up
                    //between Np and Np' if the path start to go down
                    if (shortestPath.path.get(1).getParent()==Np_pp) {
                        //going up
                        branchDistMatrix.append((shortestPath.branchDistance+Np_bl));
                    } else {
                        //going down
                        branchDistMatrix.append((shortestPath.branchDistance+Np_p_bl));
                    }
                }
            }
            nodeDistMatrix.append("\n");
            branchDistMatrix.append("\n");
            
            ////////////////////////////////////////////////////////////////////
            //fourth, build Rx virtual read datasets from removedLeaves
                 
            //Save this pruned leafs in Gx, which can be used later for 
            //simulating even more reads
            System.out.println("Write: "+GxFile.getAbsolutePath());
            FileWriter fwGx=new FileWriter(GxFile);
            for (Iterator<String> it = leavesRemoved.iterator(); it.hasNext();) {
                String next = it.next();
                fwGx.append(align.getFasta(next, false).getFormatedFasta());
                fwGx.append("\n");
            }
            fwGx.close();
                    
            for (int j = 0; j < R.length; j++) {
                System.out.println("Preparing "+(j+1)+"th query read length in R="+Arrays.toString(R)+" ; Rsd="+Rsd);
                //prepare the corresponding output files
                String readFileString=i+"_r"+R[j]+".fasta";
                File Rxj=new File(RxDir+File.separator+readFileString);
                //register it
                if (!readFiles.containsKey(i))
                    readFiles.put(i,new ArrayList<String>(R.length));
                readFiles.get(i).add(Rxj.getName());
                System.out.println("Write: "+Rxj.getAbsolutePath());
                FileWriter fwRxj=new FileWriter(Rxj);                
                int r = R[j];
                for (Iterator<String> it = leavesRemoved.iterator(); it.hasNext();) {
                    String next = it.next();
                    String seqNoGaps=align.getFasta(next, false).getSequence(true);
           
                    //build leaf_length/readLength virtual reads
                    //if virtual read length > to seqLength, use full seqLength
                    int readCount=seqNoGaps.length()/r;
                    
                    if (readCount==0) { // because readLength inferior to r
                        int p=0;
                        String read=seqNoGaps;
                        Fasta n=new Fasta(next+"_r"+R[j]+"_"+p+"_"+read.length(), read);
                        fwRxj.append(n.getFormatedFasta()+"\n"); 
                    }
                    
                    for (int k = 0; k < readCount; k++) {
                        //select normally distibuted read length v_length
                        //centered around R[j] and with standard dev of value Rsd
                        //mySample = r.nextGaussian()*desiredStandardDeviation+desiredMean
                        //The mean of the sample point is 0, and the standard deviation is 1;
                        //that means that the original sample is also its own z-score.
                        //absolute value of z represents the distance between
                        //the raw score and the population mean in units of the
                        //standard deviation.
                        //The formula is z=(x-mean)/stdev, so with the default values z=x.
                        //If we wanted to retain the z score for the sample but
                        //change the mean and stdev: 
                        //z*stdev + mean = x' where z=x, and x' represents
                        //the sample from the distribution with the desired mean
                        //and standard deviation.
                        int v_length = Double.valueOf(0.0+r+rand.nextGaussian()*Rsd).intValue();
                        //System.out.println(expString+" v_length= "+v_length);
                        
                        //if generated length is smaller than sequence length
                        if (v_length>=Rmin) {
                            //select start position 
                            int pStart= k*r;
                            //select end position, v_length or less if overcome last position
                            int pEnd= pStart+v_length;
                            if (pEnd>seqNoGaps.length()) { pEnd=seqNoGaps.length(); }
                            String read=seqNoGaps.substring(pStart, pEnd);
                            Fasta n=new Fasta(next+"_r"+R[j]+"_"+pStart+"_"+pEnd, read);
                            fwRxj.append(n.getFormatedFasta()+"\n");
                        }
                        

                    }

                    
//                    //below previous algo, where a single read per leaf was
//                    //extract randomly along leaf seqeunce
//                    
//                    int v_length = new Double(0.0+r+rand.nextGaussian()*Rsd).intValue();
//                    //if generated length is smaller than sequence length
//                    if (v_length>=Rmin && v_length<seqNoGaps.length()) {
//                        //select an uniformely distibution position in [0,seqLength-v_length]
//                        int pStart=rand.nextInt(seqNoGaps.length()-v_length+1);//this is an exclusive upper bound, so +1 to include last possibility
//                        String read=seqNoGaps.substring(pStart, pStart+v_length);
//                        Fasta n=new Fasta(next+"_r"+R[j]+"_"+pStart+"_"+(pStart+v_length-1), read);
//                        fwRxj.append(n.getFormatedFasta()+"\n");
//                    //if generated length longer than sequence, take whole sequence
//                    } else if (v_length>=Rmin) {
//                        int pStart=0;//this is an exclusive upper bound, so +1 to include last possibility
//                        String read=seqNoGaps;
//                        Fasta n=new Fasta(next+"_r"+R[j]+"_"+pStart+"_"+read.length(), read);
//                        fwRxj.append(n.getFormatedFasta()+"\n");
//                    }
                    
                    
                }
                fwRxj.close();
                
            }
            
            //know can free memory for Nx, not used anymore
            Nx.removeFromParent();
            Nx=null;
            //debug
            //tree.displayTree();

            
            
            ////////////////////////////////////////////////////////////////////
            //closes everything
            alignCopy=null;
            rootCopy=null;
            treeCopy=null;
            //do not set to null rootCopy2/treecopy2,
            //as it backs the extendedtree which will be serialized
            
            // a garbage call to clean a bit before next pruning
            System.gc();
            
            System.out.println("");
            
        } //END of nx_nodeId loop
        
        
        
        
        
        
        
        
        
        System.out.println("############################################");
        System.out.println("# Actually pruned:"+actuallyPrunedNodeIdsIndexes.size());
        System.out.println("Actually pruned, prunedNodeIds indexes:"+actuallyPrunedNodeIdsIndexes);
        System.out.print("Actually pruned,  prunedNodeIds value:");
        actuallyPrunedNodeIdsIndexes.forEach((index)->{System.out.print(" "+prunedNodeIds[index.intValue()]);});
        System.out.println("");
        System.out.println("Will be stored in expect_placement.bin such as");
        System.out.println("NxIndex map(nodeId)=index: "+NxIndex);
        System.out.println("Dtx and D'tx size (x,y): "+nodeIds.length+","+NxIndex.size());
        System.out.println("Expected placements:"+expectedPlacements);
        System.out.println("prunedTrees size : "+prunedTrees.size());
        System.out.println("############################################");
        

        
        //after all placements, save expected placements in binary file
        oos.writeObject(NxIndex);
        oos.writeObject(pruningIndex);
        oos.writeObject(expectedPlacements);
        oos.writeObject(prunedTrees);
        oos.writeObject(prunedTreesTrifurcations);
            
       
        //closes all I/O 
        brDtx.append(nodeDistMatrix);
        brD2tx.append(branchDistMatrix);
        brDtx.close();
        brD2tx.close();
        skippedLogBw.close();
        oos.close();
         
    }
    


    
    
    /**
     * shuffle an int array
     * @param array 
     */
    private void shuffleArray(Integer[] array) {
        for (int i = 0; i < array.length; i++) {
               array[i]=i;
            }
            Random generator = null;
            if (seed!=null) {
                generator=new Random(seed);
            } else {
                generator=new Random(System.nanoTime());
            }
            for (int i = 0; i < array.length - 1; i++) {
              int j = i + generator.nextInt(array.length - i);
              int t = array[j];
              array[j] = array[i];
              array[i] = t;
        }
    }
    
 
    
}
