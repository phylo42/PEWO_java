
import alignement.Alignment;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.PosixFilePermission;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import main_v2.ARProcessLauncher;
import org.jfree.io.IOUtils;
import tree.ExtendedTree;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.PhyloTree.Path;
import tree.PhyloTreeModel;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * takes Ax tree and alignment as input, then build all the directories necessary.
 * Produces also the AR commands as qsub command list
 * to the method comparisons
 * @author ben
 */
public class PrunedTreeGenerator {
    
    Long seed=new Long(1);
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");
    
    //default dataset
    File alignFile=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set.aln");
    File treeFile=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set.tree");
    
    //nodes selected for pruning at this launch
    Integer[] prunedNodeIds=null;
    //nodes effectively pruned because didn'Tx produce tre of less than 4 leaves
    ArrayList<Integer> actuallyPrunedNodeIdsIndexes=new ArrayList<>();
    
    //list of new files
    public ArrayList<File> prunedAlignmentsFiles=new ArrayList<>(); //list of pruned alignments
    public ArrayList<File> prunedTreesFiles=new ArrayList<>(); //list of pruned Trees
    public HashMap<String,List<File>> readFiles=new HashMap<>(); //map of reads, map(A0_nx4_la)=[A0_nx4_la_r150,A0_nx4_la_r300]
    public ArrayList<PhyloTree> prunedTrees=new ArrayList<>(); //list of pruned Trees
    public ArrayList<Alignment> prunedAlignments=new ArrayList<>(); //list of pruned Trees
    public File fileDtx=null;
    public File fileD2tx=null;
    
    //pruning percent
    double percentPruning=0.5; //10%

    //read generation: nomral distrib around mean R with sd (R/4)
    //and min length m
    int[] R={3,6};
    int m=2; //let's consider that we have at least 75bp reads

    //set which k/alpha are tested (1 directory created par combination
    int k=5;
    int maxK=12;
    int kIncrement=1;
    double factor=1.0;
    double maxFactor=2.0;
    double factorIncrement=0.1;
    
    //set the critera for branch injection
    float minBranchLength=-1.0f;
    int branchPerLength=1;
    
    //AR program
    File ARExecutablePath=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/paml4.9b_hacked/bin/baseml");
    //hmm programs
    File HMMBinariesDir=new File("/home/ben/Dropbox/viromeplacer/test_datasets/software/hmmer-3.1b2/binaries");
    File HMMBUILDPath=new File(HMMBinariesDir.getAbsolutePath()+File.separator+"hmmbuild");
    File HMMALIGNPath=new File(HMMBinariesDir.getAbsolutePath()+File.separator+"hmmalign");
    
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
        
        System.out.println("ARGS: workDir ARBinaries HMMBinariesDir align tree percentPruning readSize1,readSize2,...");
        
        try {
            //launch
            PrunedTreeGenerator ptg=new PrunedTreeGenerator();
            
            //PREPARE DATA SOURCES
            ///////////////////////////////////////////////////
            if (args.length>0) {
                ptg.workDir=new File(args[0]);
                ptg.ARExecutablePath=new File(args[1]);
                ptg.HMMBinariesDir=new File(args[2]);
                ptg.HMMBUILDPath=new File(ptg.HMMBinariesDir.getAbsolutePath()+File.separator+"hmmbuild");
                ptg.HMMALIGNPath=new File(ptg.HMMBinariesDir.getAbsolutePath()+File.separator+"hmmalign");
                ptg.alignFile=new File(args[3]);
                ptg.treeFile=new File(args[4]);
                ptg.percentPruning=Double.parseDouble(args[5]);
                String[] readSizes=args[6].split(",");
                ptg.R=new int[readSizes.length];
                for (int i = 0; i < readSizes.length; i++) {
                    ptg.R[i]=Integer.valueOf(readSizes[i]);
                }
            } //if no args use default values
            
            
            
            System.out.println("workDir: "+ptg.workDir);
            System.out.println("ARExecutable: "+ptg.ARExecutablePath);
            System.out.println("HMMAlignExecutable: "+ptg.HMMALIGNPath);
            System.out.println("HMMBuildExecutable: "+ptg.HMMBUILDPath);
            System.out.println("alignFile: "+ptg.alignFile);
            System.out.println("treeFile: "+ptg.treeFile);
            System.out.println("readSizes: "+Arrays.toString(ptg.R));

            //TEST ZONE
            
            //very small default dataset for debugging
//            File treeFile=new File(ptg.workDir.getParent()+File.separator+"6_leaves_test_set.tree");
//            File alignFile=new File(ptg.workDir.getParent()+File.separator+"6_leaves_test_set.aln");
//            ptg.percentPruning=1.0;
            //pplacer rRNA dataset
            
            //directory where all testes are done
//            File workDir=new File("/media/ben/STOCK/DATA/viromeplacer/accu_tests");
//            File treeFile=new File(workDir+File.separator+"RAxML_result.bv_refs_aln");
//            File alignFile=new File(workDir+File.separator+"bv_refs_aln_stripped_99.5.fasta");
//            //set directory associated to this tree
//            ptg.workDir=new File(workDir.getAbsolutePath()+File.separator+"pplacer_16s");


//            
//            //read tree
//            BufferedReader br=new BufferedReader(new FileReader(treeFile));
//            String line=null;
//            String treeString=null;
//            while((line=br.readLine())!=null) {treeString=line;}
//            PhyloTree tree = NewickReader.parseNewickTree2(treeString, true);
//            tree.initIndexes();
//            System.out.println("Original tree, # nodes: "+tree.getNodeCount());
//            //tree.displayTree();

                    
            //LOAD TREE / ALIGNMENTS
            ///////////////////////////////////////////////////
            //load alignment
            FASTAPointer fp=new FASTAPointer(ptg.alignFile, false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
            }
            Alignment align=new Alignment(fastas);
            System.out.println(align.describeAlignment(false));
            fp.closePointer();
            //load tree
            BufferedReader br=new BufferedReader(new FileReader(ptg.treeFile));
            String line=null;
            String treeString=null;
            while((line=br.readLine())!=null) {treeString=line;}
            PhyloTree tree = NewickReader.parseNewickTree2(treeString, true, false);
            tree.initIndexes();
            System.out.println("Original tree, # nodes: "+tree.getNodeCount());
            //tree.displayTree();
            
            
            
            
            //PREPARE ALL PRUNING EXPERIMENT FILES (Ax, Tx, Rx)
            //Ax: the pruned alignments
            //Tx: the pruned trees
            //Rx: the reads built from the pruned leaves
            ///////////////////////////////////////////////////
            if (!ptg.workDir.exists()) {
                ptg.workDir.mkdir();
            }
            ptg.generatePrunedTrees(ptg.workDir,tree,align);
            System.out.println("PRUNED TREE GENERATION DONE !");
            
            
            
            
            
            //PREPARE THE AR direcoties FOR ALL k/alpha COMBINATIONS (Dx)
            //i.e : build the extended trees
            //      build AR command, without execution
            //      prepare script which launch these using SGE (qsub)  
            ///////////////////////////////////////////////////
            ptg.prepareARcommands(ptg.workDir);
            System.out.println("BUILD OF ALL DB (different k/alpha) DONE !");
            
            
            ////////////////////////////////////////////////////////////////////
            //build HMM alignments for this current pruning experiment
            //i.e. : build the hmm profile from the corresponding Ax
            //       align all Rx reads to this profile
            //       convert the result from stockolm to fasta alignment
            ptg.prepareHMMCommands(ptg.workDir);
            System.out.println("BUILD OF ALL HMM-BASED ALIGNMENTS DONE !");
            
            
            
            Thread.sleep(1);
            
        } catch (IOException ex) {
            Logger.getLogger(PrunedTreeGenerator.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(PrunedTreeGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    
    /**
     * main script
     * @param treeCopy 
     */
    private void generatePrunedTrees(File workDir,PhyloTree tree,Alignment align) throws IOException {
        

        
        //get list of internal Nx ids
        Integer[] nodeIds=new Integer[tree.getNodeIdsByDFS().size()];
        System.out.println("This tree contains "+nodeIds.length+" nodeIds.");
        //shuffle their order
        shuffleArray(nodeIds);
        //define first x% as pruning experiments
        prunedNodeIds=Arrays.copyOfRange(nodeIds, 0, new Double(percentPruning*nodeIds.length).intValue());
        System.out.println("# pruning attempts: "+prunedNodeIds.length);
        System.out.println("prunedNodeIds, selected for pruning attempts (raw): "+Arrays.toString(prunedNodeIds));
        //sort to make more comprehensive output matrices Dtx and D'tx
        Arrays.sort(prunedNodeIds);
        System.out.println("prunedNodeIds, selected for pruning attempts (ordered): "+Arrays.toString(prunedNodeIds));

        //log registering all skipped Nx, and the reason
        File skippedLog=new File(workDir+File.separator+"SKIPPED_Nx");
        BufferedWriter skippedLogBw=new BufferedWriter(new FileWriter(skippedLog,false));
        
        //write in a binary file the pruned tree and the expected placement,
        //that is the branch b_new (and b_new_p if rerooting), see algo below.
        //expected placement is saved through an integer array of 1 or 2 elements
        //integer is the n of the node son of b_new
        File expectedPlacementsFile=new File(workDir.getAbsolutePath()+File.separator+"expected_placements.bin");
        ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(expectedPlacementsFile),4096));
        HashMap<Integer,Integer> NxIndex=new HashMap<>(); //map(120,i.e. nx120)=0; map(1142)=1 ; map(5454)=2 ...
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
        
        
            
        //launch pruning for each selected Nx
        for (int i = 0; i < prunedNodeIds.length; i++) {
            Integer nx_id = prunedNodeIds[i];
            //System.out.println("--------------------------------------");
            //System.out.println("copying tree and alignment before pruning...");
            PhyloNode rootCopy=tree.getRoot().copy();
            PhyloTree treeCopy=new PhyloTree(new PhyloTreeModel(rootCopy),tree.isRooted(), false);;
            //System.out.println("indexing tree ...");
            treeCopy.initIndexes();
//            //some checkup about the original  copy
//            System.out.println("Tree; #nodes="+treeCopy.getNodeCount());
//            System.out.println("Tree; rooted="+treeCopy.isRooted());
//            System.out.println("Tree; #leaves="+treeCopy.getLeavesCount());
//            Enumeration depthFirstEnumeration = treeCopy.getRoot().depthFirstEnumeration();
//            while(depthFirstEnumeration.hasMoreElements()) {
//                System.out.println("Tree;  nodes="+depthFirstEnumeration.nextElement());
//            }
            //copying alignment
            Alignment alignCopy=align.copy();
            //System.out.println("Starting pruning...");
            
            //current root Nx defining the pruned clade
            PhyloNode Nx= treeCopy.getById(nx_id);
            System.out.println("--------------------------------------");
            System.out.println("selected Nx: x="+i+"  node:"+Nx);
            if (Nx.isRoot()) {
                skippedLogBw.append("Nx="+i+"\t"+Nx.toString()+"\tIs root, so skipped.\n");
                System.out.println("SKIPPED: this node is root (id="+nx_id+"), no pruning.");
                continue;
            }
                        
            
            
            
            ///////////////////////////////////////////////////////////
            //first, modify alignment and deleted subtree related to this Nx
            
            
            //if leaf, removes from multiple align
            //and export this leave to Ax it file
            ArrayList<Fasta> leavesRemoved=new ArrayList<>();
            if (Nx.isLeaf()) {
                leavesRemoved.add(alignCopy.getFasta(Nx.getLabel(), false));
                alignCopy.removeSequence(Nx.getLabel());
            //if an internal Nx
            } else {
                //enumerate nodes in subtree
                PhyloNode nextNx =null;
                Enumeration<PhyloNode> DFSenum = Nx.depthFirstEnumeration();
                //careful, tree topology change are reflected in the enumeration
                //so we need to remove node AFTER this this transversal postorder
                while (DFSenum.hasMoreElements()) {
                    nextNx=DFSenum.nextElement();
                    //if leaf, removes and export
                    if (nextNx.isLeaf()) {
                        leavesRemoved.add(alignCopy.getFasta(nextNx.getLabel(), false));
                        alignCopy.removeSequence(nextNx.getLabel());
                    } 
                }
                //delete Nx children (this does not delete the sub nodes from memory !!!)
                nextNx.removeAllChildren();
            }
            //if this removal ends to Ax multiple alignment with only <3 leaves
            //we don'Tx go further and pass to the next Nx.
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
            
            
            //if reach here, we made Ax tree with more than 3 leaves, this pruning
            // will be effectively operated
            actuallyPrunedNodeIdsIndexes.add(i);
            
            //so let's register this pruning index in the expected_placement.bin
            //file
            NxIndex.put(nx_id,actuallyPrunedNodeIdsIndexes.size()-1);
            
            
            //prepare the corresponding output files
            File AxDir=new File(workDir+File.separator+"Ax");
            File TxDir=new File(workDir+File.separator+"Tx");
            File RxDir=new File(workDir+File.separator+"Rx");
            AxDir.mkdir();
            TxDir.mkdir();
            RxDir.mkdir();
            
            File Ax=new File(AxDir+File.separator+"A"+i+"_nx"+nx_id+"_la"+Nx.getLabel()+".align");
            File Tx=new File(TxDir+File.separator+"T"+i+"_nx"+nx_id+"_la"+Nx.getLabel()+".tree");
            prunedTreesFiles.add(Tx);
            prunedAlignmentsFiles.add(Ax);
            
            System.out.println("LeavesRemoved: "+leavesRemoved);
            
            ////////////////////////////////////////////////////////////////////
            //second, built Tx.
            //removal of node Nx of Tx is operated on the basis of the following
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

            //memorize which edge is the best placement, this is Ax node, 
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
                treeCopy=new PhyloTree(new PhyloTreeModel(Np),tree.isRooted(), false);
                //memorize best placement
                b_new=Np_p;
                b_new_p=Np_pp;
                

                
            //in all other case, the tree is disconnected and reconnected
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
            alignCopy.writeAlignmentAsFasta(Ax);
            prunedAlignments.add(alignCopy);
            //save the pruned tree
            System.out.println("Indexing pruned tree");
            treeCopy.initIndexes(); //necessary to remove former nodes from the maps shortcuts
            System.out.println("Writing pruned tree newick");
            NewickWriter nw=new NewickWriter(Tx);            
            nw.writeNewickTree(treeCopy, true, true, false);  //no internal node names if PAML !
            prunedTrees.add(treeCopy);
            
            //SAVE THE EXPECTED PLACEMENT IN A FILE,
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
//            System.out.println(Arrays.toString(prunedNodeIds));
//            System.out.println("Np_p:"+Np_p);
//            System.out.println("Np_pp:"+Np_pp);
//            System.out.println("b_new:"+b_new);
            //using the nodeIds table, Dtx column order will match the
            //shuffled nodeIds, like this the first xx% correspond to the
            //pruned nodes. we could have used tree.getNodeIdsByDFS() too.
            for (int n=0;n<nodeIds.length;n++) {
                PhyloNode currentNode=treeCopy.getById(nodeIds[n]);
                //System.out.println("+++++++++++currentNode:"+currentNode);
                nodeDistMatrix.append(";");
                branchDistMatrix.append(";");
                if (currentNode==null) { //this node was pruned, so absent from treeCopy
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
                    //added_root comes from forced rooting of the input tree
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
            
          
//            System.out.println("----------------");
//            System.out.println(nodeDistMatrix.toString());
//            System.out.println("----------------");
//            System.out.println(branchDistMatrix.toString());
//            System.out.println("----------------");

            ////////////////////////////////////////////////////////////////////
            //fourth, build virtual read datasets from removedLeaves

            String expString="R"+i+"_nx"+nx_id+"_la"+Nx.getLabel();
            for (int j = 0; j < R.length; j++) {
                //prepare the corresponding output files
                String readFileString=expString+"_r"+R[j]+".fasta";
                String experimentAlignmentLabel=Ax.getName().split("\\.")[0];
                File Rxj=new File(RxDir+File.separator+readFileString);
                //register it
                if (!readFiles.containsKey(experimentAlignmentLabel))
                    readFiles.put(experimentAlignmentLabel,new ArrayList<File>(R.length));
                readFiles.get(experimentAlignmentLabel).add(Rxj);
                       
                FileWriter fw=new FileWriter(Rxj);
                int r = R[j];
                for (Iterator<Fasta> it = leavesRemoved.iterator(); it.hasNext();) {
                    Fasta next = it.next();
                    String seqNoGaps=next.getSequence(true);
                    //select normally distibuted read length v_length
                    //centered around R[j] and with standard dev of r/4
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
                    Random rand = new Random(seed);
                    int v_length = new Double(r+rand.nextGaussian()*(new Double(r)/4.0)).intValue();
                    //System.out.println(v_length);
                    if (v_length>=m && v_length<seqNoGaps.length()) {
                        //select an uniformely distibution position in [0,seqLength-v_length]
                        int p=rand.nextInt(seqNoGaps.length()-v_length+1);//this is an exclusive upper bound, so +1 to includelast possibility
                        String read=next.getSequence(true).substring(p, p+v_length);
                        Fasta n=new Fasta(next.getHeader()+"_r"+R[j]+"_"+p+"_"+(p+v_length-1), read);
                        fw.append(n.getFormatedFasta()+"\n");
                    }
                }
                fw.close();
                
            }
            
            //know can free memory for Nx, not used anymore
            Nx.removeFromParent();
            Nx=null;
            //debug
            //tree.displayTree();
            
            //closes everything
            nw.close();
            treeCopy=null;
                        
            
            
        }
        
        
        System.out.println("############################################");
        System.out.println("# Actually pruned:"+actuallyPrunedNodeIdsIndexes.size());
        System.out.println("Actually pruned (prunedNodeIds indexes):"+actuallyPrunedNodeIdsIndexes);
        System.out.println("Will be stared in expect_placement.bin such as");
        System.out.println("NxIndex map(nodeId)=index: "+NxIndex);
        System.out.println("Dtx and D'tx size (x,y): "+nodeIds.length+","+NxIndex.size());
        System.out.println("############################################");
        
        
        //after all placements, save expected placements in binary file
        oos.writeObject(NxIndex);
        oos.writeObject(expectedPlacements);
        oos.writeObject(prunedTrees);
        
        
        
        
        //closes all I/O 
        brDtx.append(nodeDistMatrix);
        brD2tx.append(branchDistMatrix);
        brDtx.close();
        brD2tx.close();
        skippedLogBw.close();
        oos.close();
         
    }
    
    
    /**
     * main script
     * @param workDir 
     * @throws IOException 
     */
    private void prepareARcommands(File workDir) throws IOException {

        //to ensure consistent direcory names (alpha-> _ax_)
        NumberFormat nf = NumberFormat.getNumberInstance();
        nf.setMinimumFractionDigits(1);
        nf.setMaximumFractionDigits(1);
        //
        ArrayList<File> workDirs=new ArrayList<>();
        
        
        //all AR reconstructions and viromeplacer DBs in Dx
        File Dx=new File(workDir.getAbsolutePath()+File.separator+"Dx");
        Dx.mkdir();
        //File where to save all qsub commands
        File qSubCommands=new File(workDir+File.separator+"Dx"+File.separator+"qsub_AR_commands");
        BufferedWriter bw=new BufferedWriter(new FileWriter(qSubCommands));
        
        
        
        
        //launch pruning for each selected Nx
        for (int i = 0; i < actuallyPrunedNodeIdsIndexes.size(); i++) {
            Integer nx_id = actuallyPrunedNodeIdsIndexes.get(i);
            File a=prunedAlignmentsFiles.get(i);
            File t=prunedTreesFiles.get(i);
            
            String experimentAlignmentLabel=a.getName().split("\\.")[0];
            String experimentTreeLabel=t.getName().split("\\.")[0];
            
            System.out.println("PREPARING ALL DB_BUILDs FOR EXPERIMENT: "+experimentAlignmentLabel);

            //alignment filename is used as default workDir for this experiment
            File expPath=new File(Dx+File.separator+experimentAlignmentLabel);
            System.out.println("Experiment work dir: "+expPath.getAbsolutePath());
            workDirs.add(expPath);
            if (workDir.canWrite())
                expPath.mkdir();
            else {
                System.out.println("Cannot write in dir: "+expPath.getAbsolutePath());
                System.exit(1);
            }
            
            
            //load data necessary to do the extended tree
            ///////////////////////////////////////////////////
            
            //directory of extended tree
            File extendedTreeDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"extended_trees");
            extendedTreeDir.mkdir();
            //select the tree and alignments of current pruning experiment
            // workDir/Tx/TX_nxX_laX.tree
            File Tx=new File(workDir.getAbsolutePath()+File.separator+"Tx"+File.separator+experimentTreeLabel+".tree");
            // workDir/Ax/AX_nxX_laX.align
            File Ax=new File(workDir.getAbsolutePath()+File.separator+"Ax"+File.separator+experimentAlignmentLabel+".align");
            //load alignment of this experiment 
            Alignment experimentAlign=prunedAlignments.get(i);
            System.out.println(experimentAlign.describeAlignment(false));
            //load tree of this experiment
            PhyloTree experimentTree = prunedTrees.get(i); //indexes already init at creation
            System.out.println("Experiment tree, # nodes: "+experimentTree.getNodeCount());            
            //build the extended trees first and save them in the extended_trees directory
            //these will be used in the AR below
            System.out.println("Loading extended tree from binary version.");
            ExtendedTree experimentExtendedTrees = buildExtendedTrees(extendedTreeDir, experimentAlign, experimentTree);
            
            
            
            File fileRelaxedAlignmentPhylip=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.phylip");
            File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
            File fileRelaxedAlignmentFasta=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.fasta");
            File fileRelaxedTreewithBL=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
            File fileRelaxedTreeBinary=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree.bin");
            
            //prepare the directory for AR based on extended tree 
            ///////////////////////////////////////////////////
                    
            //create correpsonding directory        
            File ARDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"AR");
            ARDir.mkdir();
            System.out.println("   AR in: "+ARDir.getAbsolutePath());
            ARProcessLauncher arpl=new ARProcessLauncher(ARProcessLauncher.AR_PAML, ARExecutablePath, false);
            
            //prepare corresponding AR commands (with PAML)
            String ARCommand = arpl.prepareAR(ARDir, fileRelaxedAlignmentPhylip, fileRelaxedTreewithBLNoInternalNodeLabels);
            bw.append("echo \""+ARCommand+"\" | qsub -N "+expPath.getName()+" -wd "+ARDir.getAbsolutePath());
            bw.newLine();

            //make subdirectory corresponding to  k/alpha combinations
            //this will be the directory in which viromeplacer will work
            for (int current_k = k; current_k < maxK+kIncrement; current_k+=kIncrement) {
                for (double current_alpha = factor; current_alpha < maxFactor+factorIncrement; current_alpha+=factorIncrement) {

                    Infos.println("#######################################################");
                    Infos.println("# PREPARING COMMANDS FOR QSUB AR LAUNCHES; Config: k="+current_k+" factor="+current_alpha);
                    Infos.println("#######################################################");

                    //make subdirectory corresponding to this k/alpha combination
                    File combDir=new File(expPath.getAbsolutePath()+File.separator+"k"+current_k+"_a"+nf.format(current_alpha));
                    boolean dirOK=combDir.mkdir();
                    //do Ax subdir AR/ in each k/alpha directory
                    //make symbolic link to rst file which is in experiment directory (same for all k/alpha combinations)
                    File combDirAR=new File(combDir.getAbsolutePath()+File.separator+"AR");
                    combDirAR.mkdir();
                    File combDirExtendedTree=new File(combDir.getAbsolutePath()+File.separator+"extended_trees");
                    combDirExtendedTree.mkdir();
                    File originalRST=new File(ARDir.getAbsolutePath()+File.separator+"rst");
                    originalRST.createNewFile();
                    File linkToRST=new File(combDirAR.getAbsolutePath()+File.separator+"rst");
                    if (!linkToRST.exists())
                        Files.createSymbolicLink(linkToRST.toPath(), originalRST.toPath());
                    //do Ax subdir extended_trees/ in each k/alpha directory
                    //make symbolic link to correpsonding trees/alignments
                    //which are in experiment directory (same trees/aligns for all k/alpha combinations)
                    //this is necessary because these trees are required by db_build
                    File linkToTree=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
                    File linkToTreeNoInter=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
                    File linkToAlignmentFasta=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_align.fasta");
                    File linkToAlignmentPhylip=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_align.phylip");
                    File linkToRelaxedTreeBinary=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree.bin");
                    
                    if (!linkToTree.exists())
                        Files.createSymbolicLink(linkToTree.toPath(), fileRelaxedTreewithBL.toPath());       
                    if (!linkToTreeNoInter.exists())
                        Files.createSymbolicLink(linkToTreeNoInter.toPath(), fileRelaxedTreewithBLNoInternalNodeLabels.toPath());  
                    if (!linkToAlignmentFasta.exists())
                        Files.createSymbolicLink(linkToAlignmentFasta.toPath(), fileRelaxedAlignmentFasta.toPath()); 
                    if (!linkToAlignmentPhylip.exists())
                        Files.createSymbolicLink(linkToAlignmentPhylip.toPath(), fileRelaxedAlignmentPhylip.toPath());
                    if (!linkToRelaxedTreeBinary.exists()) 
                        Files.createSymbolicLink(linkToRelaxedTreeBinary.toPath(), fileRelaxedTreeBinary.toPath());
                    

                    
                }

            }



        }
        
        //give execution rights to the commands file
        Files.setPosixFilePermissions(qSubCommands.toPath(), perms);
        bw.close();
        
        
    }
    
    
    /**
     * main script
     * @param workDir 
     * @throws IOException 
     */
    private void prepareHMMCommands(File workDir) throws IOException {
        
        //prepare files with all commands necessary for the HMM alignment build
        //1. commands for 1 alignment are saved as Ax "set" of commands
        //in Ax sh script.
        //2. all sh script calls are associated to qsub commands and put in
        //a single file
        
        StringBuilder qSubCommands=new StringBuilder();
            
        //create HMMx directories
        File HMMxDir=new File(workDir+File.separator+"HMMx");
        HMMxDir.mkdir();
        //writee a copy of the psiblast2fasta.py script from the jar to the HMMxDir
        File pyScript=new File(HMMxDir.getAbsolutePath()+File.separator+"psiblast2fasta.py");
        InputStream resourceAsStream = this.getClass().getClassLoader().getResourceAsStream("scripts/psiblast2fasta.py");
        Files.copy(resourceAsStream, pyScript.toPath(), StandardCopyOption.REPLACE_EXISTING);
        resourceAsStream.close();
        Files.setPosixFilePermissions(pyScript.toPath(), perms);

        
        //launch pruning for each selected Nx
        for (int i = 0; i < actuallyPrunedNodeIdsIndexes.size(); i++) {
            Integer nx_id = actuallyPrunedNodeIdsIndexes.get(i);
            File Ax=prunedAlignmentsFiles.get(i);
            File Tx=prunedTreesFiles.get(i);
            String experimentAlignmentLabel=Ax.getName().split("\\.")[0];
            String experimentTreeLabel=Tx.getName().split("\\.")[0];
            System.out.println("PREPARING HMM ALIGNEMENTS FOR EXPERIMENT: "+experimentAlignmentLabel);

            //alignment filename is used as default workDir for this experiment
            File HMMxAxDir=new File(HMMxDir.getAbsolutePath()+File.separator+experimentAlignmentLabel);
            System.out.println("Experiment work dir: "+HMMxAxDir.getAbsolutePath());
            if (workDir.canWrite())
                HMMxAxDir.mkdir();
            else {
                System.out.println("Cannot write in dir: "+HMMxAxDir.getAbsolutePath());
                System.exit(1);
            }
            
            //all files used in this HMM alignment
            File alignment=Ax;
            File hmmOuput=new File(HMMxAxDir.getAbsolutePath()+File.separator+experimentAlignmentLabel+".hmm");
            
            //build the hmmbuild command. 1 hmm per experiment
            StringBuilder commandSet=new StringBuilder();
            //commandSet.append("#!/bin/sh\n");
            commandSet.append(HMMBUILDPath.getAbsolutePath()+" --dna "+hmmOuput.getAbsolutePath()+" "+Ax.getAbsolutePath()+"\n");
            
            //build alignments of reads
            List<File> reads=readFiles.get(experimentAlignmentLabel);
            for (int j=0;j<reads.size();j++) {
                File alnOutput=new File(HMMxAxDir.getAbsolutePath()+File.separator+reads.get(j).getName().split("\\.")[0]+".aln.psiblast");
                commandSet.append(HMMALIGNPath.getAbsolutePath()+" --outformat PSIBLAST "
                                    + " --mapali "+Ax.getAbsolutePath()+" -o "+alnOutput.getAbsolutePath()
                                    + " "+hmmOuput+" "+reads.get(j)
                                    + "\n");     
                //conversion to fasta
                commandSet.append(  pyScript.getAbsolutePath()+" "+
                                    alnOutput.getAbsolutePath()+" "+
                                    alnOutput.getParentFile().getAbsolutePath()+File.separator+alnOutput.getName().split("\\.")[0]+".aln.fasta" +
                                    "\n");   
            }
            //save this alignment commands in a sh script
            File shScript=new File(HMMxAxDir.getAbsolutePath()+File.separator+"align_by_hmm.sh");
            BufferedWriter bwShScript=new BufferedWriter(new FileWriter(shScript));
            bwShScript.append(commandSet.toString());
            bwShScript.close();  
            Files.setPosixFilePermissions(shScript.toPath(), perms);
            //add the qsub execution line of this sh script in the list of qsub
            //commands
            qSubCommands.append("echo \""+shScript.getAbsolutePath()+"\" |  qsub -N HMMx_"+experimentAlignmentLabel+
                                " -wd "+HMMxAxDir.getAbsolutePath()+"\n");
        }
        //write the qsub command file
        File qsubScript=new File(HMMxDir.getAbsolutePath()+File.separator+"qsub_hmm_commands");
        BufferedWriter bwQsubScript=new BufferedWriter(new FileWriter(qsubScript));
        bwQsubScript.append(qSubCommands.toString());
        bwQsubScript.close();
        Files.setPosixFilePermissions(qsubScript.toPath(), perms);
        

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
    
    /**
     * execution itself
     * @param com 
     */
    private static void executeProcess(List<String> com, File workDir) throws IOException {

        ProcessBuilder pb = new ProcessBuilder(com);
        //pb.environment().entrySet().stream().forEach((e) ->{ System.out.println(e.getKey()+"="+e.getValue()); });
        //env.put("VAR1", "myValue"); env.remove("OTHERVAR");
        pb.directory(workDir);
        Infos.println("Current directory:"+pb.directory().getAbsolutePath());
        pb.redirectErrorStream(false);
        pb.redirectOutput(ProcessBuilder.Redirect.PIPE);
        pb.redirectInput(ProcessBuilder.Redirect.PIPE);
        Process p = pb.start();
        assert pb.redirectInput() == ProcessBuilder.Redirect.PIPE;
        assert p.getInputStream().read() == -1;
        //redirect sdtout/stdin to files
        FileOutputStream STDOUTOutputStream=new FileOutputStream(new File(workDir.getAbsolutePath()+File.separator+"_sdtout.txt"));
        FileOutputStream STDERROutputStream=new FileOutputStream(new File(workDir.getAbsolutePath()+File.separator+"_sdterr.txt"));
        inputStreamToOutputStream(p.getInputStream(), System.out);
        inputStreamToOutputStream(p.getInputStream(), STDOUTOutputStream);
        inputStreamToOutputStream(p.getErrorStream(), STDERROutputStream);
        Infos.println("External process operating reconstruction is logged in: "+new File(workDir.getAbsolutePath()+File.separator+"sdtout.txt").getAbsolutePath());
        Infos.println("Launching process ...");
        try {
            p.waitFor();
            Thread.sleep(1000);
        } catch (InterruptedException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.flush();
        STDOUTOutputStream.flush();
        STDERROutputStream.flush();
        STDOUTOutputStream.close();
        STDERROutputStream.close();
        System.out.println(""); //this line ensures line return after the external process output
        System.out.println("Process finished.");
    }
    
    
    public static void inputStreamToOutputStream(final InputStream inputStream, final OutputStream out) {
        Thread t = new Thread(new Runnable() {
            @Override
            public void run() {
                try {
                    int d;
                    BufferedInputStream bis=new BufferedInputStream(inputStream);
                    while ((d = bis.read()) != -1) {
                        out.write(d);
                    }
                } catch (IOException ex) {
                    Infos.println(ex.getCause());
                }
            }
        });
        t.setDaemon(false);
        t.start();
    }
    
    
    private ExtendedTree buildExtendedTrees(File extendedTreePath,Alignment align, PhyloTree tree) {
        
        ExtendedTree extendedTree = null;

        File fileRelaxedAlignmentFasta=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_align.fasta");
        File fileRelaxedAlignmentPhylip=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_align.phylip");
        File fileRelaxedTreewithBL=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
        File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");;
        File fileRelaxedTreeBinary=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_tree.bin");
        
        try {
            System.out.println("Injecting fake nodes...");
            //note, we read again the tree to build Ax new PhyloTree object
            //this is necessary as its TreeModel is directly modified
            //at instanciation of ExtendedTree
            extendedTree=new ExtendedTree(tree,minBranchLength,branchPerLength);                    
            extendedTree.initIndexes(); // don'Tx forget to reinit indexes !!!
            ArrayList<PhyloNode> listOfNewFakeLeaves = extendedTree.getFakeLeaves();
            System.out.println("RelaxedTree contains "+extendedTree.getLeavesCount()+ " leaves");
            System.out.println("RelaxedTree contains "+extendedTree.getFakeLeaves().size()+ " FAKE_X new leaves");
            //add new leaves to alignment
            char[] gapSeq=new char[align.getLength()];
            Arrays.fill(gapSeq, '-');
            ArrayList<char[]> seqs=new ArrayList<>();
            String[] labels=new String[listOfNewFakeLeaves.size()];
            for (int i = 0; i < listOfNewFakeLeaves.size(); i++) {
                labels[i]=listOfNewFakeLeaves.get(i).getLabel();
                seqs.add(gapSeq);
            }
            align.addAllSequences(labels,seqs);
            
            //write alignment and ARTree for BrB
            System.out.println("Write extended alignment (fasta): "+fileRelaxedAlignmentFasta.getAbsolutePath());
            align.writeAlignmentAsFasta(fileRelaxedAlignmentFasta);
            System.out.println("Write extended alignment (phylip): "+fileRelaxedAlignmentPhylip.getAbsolutePath());
            align.writeAlignmentAsPhylip(fileRelaxedAlignmentPhylip);
            align=null;
            //write extended trees
            System.out.println("Write extended newick tree: "+fileRelaxedTreewithBL.getAbsolutePath());
            NewickWriter nw=new NewickWriter(fileRelaxedTreewithBL);
            nw.writeNewickTree(extendedTree, true, true, false);
            nw.close();
            //write version without internal nodes labels
            System.out.println("Write extended newick tree with branch length: "+fileRelaxedTreewithBLNoInternalNodeLabels.getAbsolutePath());
            nw=new NewickWriter(fileRelaxedTreewithBLNoInternalNodeLabels);
            nw.writeNewickTree(extendedTree, true, false, false);
            //save this extendedTree as Ax binary
            FileOutputStream fos = new FileOutputStream(fileRelaxedTreeBinary);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
            Infos.println("Storing binary version of Extended Tree.");
            oos.writeObject(extendedTree);
            oos.close();
            //extendedTreeForJplace=nw.getNewickTree(extendedTree, true, true, true);
            nw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
            System.out.println("Error raised from extended tree reconstruciton!");
        }
        
        return extendedTree;
    }
    
}
