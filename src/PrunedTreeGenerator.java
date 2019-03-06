
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
import java.io.InputStream;
import java.io.ObjectOutputStream;
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
import java.util.Locale;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
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
 * From an alignment and tree, build directories of pruning experiments, 
 if tree is unrooted then roots it for each pruning experiment,
 saves pruned alignments in AxFile, pruned trees in TxFile and generated reads in Rx,
 build Dtx matrice which defines node_dist between edges (note that roots 
 are not considered in this distance, the 2 edges sons of root are considered
 as the same placement) 
 * @author ben
 */
public class PrunedTreeGenerator {
    
    //pseudorandom generator
    Long seed=new Long(1);
    Random rand = new Random(seed);
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");
    
    //default dataset
    File alignFile=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set.aln");
    File treeFile=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set.tree");
    
    //States
    States s=new DNAStatesShifted();
    
    //nodes selected for pruning at this launch
    Integer[] prunedNodeIds=null;
    //nodes effectively pruned because didn'TxFile produce tre of less than 4 leaves
    ArrayList<Integer> actuallyPrunedNodeIdsIndexes=new ArrayList<>();
    
    //list of new files
    public HashMap<String,List<String>> readFiles = new HashMap<>(); //map of reads, map(A0_nx4_la)=[A0_nx4_la_r150,A0_nx4_la_r300]
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
    
    //Generate all possible trifurcations for each pruning
    //build as many databases and placement commands.
    //The same parameters and alignment file is used.
    //trifurcations follow naming convention :
    //Tx_trifuXX_nxXX_laXX
    boolean trifurcations=false;
    ArrayList<Integer> trifurcationsNxIndexes=null;
    

    //set which minK/alpha are tested (1 directory created par combination
    int minK=5;
    int maxK=12;
    int kIncrement=1;
    float minFactor=1.0f;
    float maxFactor=2.0f;
    float factorIncrement=0.1f;
    
    //set the critera for branch injection
    float minBranchLength=-1.0f; //basically all branches, of all length (even 0)
    int branchPerEdge=1;

    //hmm programs
    File HMMBinariesDir=new File("/home/ben/Dropbox/viromeplacer/test_datasets/software/hmmer-3.1b2/binaries");
    File HMMBUILDPath=new File(HMMBinariesDir.getAbsolutePath()+File.separator+"hmmbuild");
    File HMMALIGNPath=new File(HMMBinariesDir.getAbsolutePath()+File.separator+"hmmalign");
    //RAXML for tree optimisation
    File RAXMLExecutable=new File("/home/ben/Dropbox/viromeplacer/test_datasets/software/RAxML-8.2.9/raxmlHPC-PTHREADS-SSE3");
    
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
        
        System.out.println("ARGS: workDir HMMBinariesDir RAXMLBinary align tree percentPruning(float) readSize1(int),readSize2(int),... readSD(int) [ branchPerEdge[int] kmin[int] kmax[int] kstep[int] amin[float] amax[float] astep[float] [nucl=0|prot=1] [trifurcations:-1=no|1,45,48=list of Nx to test] ] ");
        
        System.out.println("Command: "+Arrays.toString(args).replaceAll(",", " "));

        
        try {
            //launch
            PrunedTreeGenerator ptg=new PrunedTreeGenerator();
            
            //PREPARE DATA SOURCES
            ///////////////////////////////////////////////////
            if (args.length>0) {
                ptg.workDir=new File(args[0]);
                ptg.HMMBinariesDir=new File(args[1]);
                ptg.HMMBUILDPath=new File(ptg.HMMBinariesDir.getAbsolutePath()+File.separator+"hmmbuild");
                ptg.HMMALIGNPath=new File(ptg.HMMBinariesDir.getAbsolutePath()+File.separator+"hmmalign");
                ptg.RAXMLExecutable=new File(args[2]);
                ptg.alignFile=new File(args[3]);
                ptg.treeFile=new File(args[4]);
                //ptg.percentPruning=Double.parseDouble(args[5]);
                ptg.pruningCount=Integer.parseInt(args[5]);
                String[] readSizes=args[6].split(",");
                ptg.R=new int[readSizes.length];
                for (int i = 0; i < readSizes.length; i++) {
                    ptg.R[i]=Integer.valueOf(readSizes[i]);
                }
                ptg.Rsd=Integer.parseInt(args[7]);
                //optionnals, if no args use default values
                if (args.length>8) {
                    ptg.branchPerEdge=Integer.parseInt(args[8]);
                }
                if (args.length>9) {
                    ptg.minK=Integer.parseInt(args[9]);
                    ptg.maxK=Integer.parseInt(args[10]);
                    ptg.kIncrement=Integer.parseInt(args[11]);
                    ptg.minFactor=Float.parseFloat(args[12]);
                    ptg.maxFactor=Float.parseFloat(args[13]);
                    ptg.factorIncrement=Float.parseFloat(args[14]);
                    int protein=Integer.parseInt(args[15]);
                    ptg.proteinAnalysis=(protein>0);
                    if (!args[16].equals("-1")) {
                        ptg.trifurcations=true;
                        String[] trifuIndexes=args[16].split(",");
                        ptg.trifurcationsNxIndexes=new ArrayList<>(trifuIndexes.length);
                        for (String trifuIndexe : trifuIndexes) {
                            ptg.trifurcationsNxIndexes.add(Integer.valueOf(trifuIndexe));
                        }
                        
                    }
                }
                
            }
            
            
            
            System.out.println("workDir: "+ptg.workDir);
            System.out.println("RAXML Executable: "+ptg.RAXMLExecutable);
            System.out.println("HMMAlignExecutable: "+ptg.HMMALIGNPath);
            System.out.println("HMMBuildExecutable: "+ptg.HMMBUILDPath);
            System.out.println("alignFile: "+ptg.alignFile);
            System.out.println("treeFile: "+ptg.treeFile);
            //System.out.println("percentPruning: "+ptg.percentPruning);
            System.out.println("pruningCount: "+ptg.pruningCount);
            System.out.println("readSizes: "+Arrays.toString(ptg.R));
            System.out.println("branchPerEdge: "+ptg.branchPerEdge);
            System.out.println("mink:"+ptg.minK);
            System.out.println("maxk:"+ptg.maxK);
            System.out.println("incrementk:"+ptg.kIncrement);
            System.out.println("minalpha:"+ptg.minFactor);
            System.out.println("maxalpha:"+ptg.maxFactor);
            System.out.println("incrementalpha:"+ptg.factorIncrement);
            System.out.println("protein:"+ptg.proteinAnalysis);
            System.out.println("trifurcations:"+ptg.trifurcations);
            System.out.println("trifurcations Nx tested: "+ptg.trifurcationsNxIndexes);

            //TEST ZONE
            
            //very small default dataset for debugging
//            File treeFile=new File(ptg.workDir.getParent()+File.separator+"6_leaves_test_set.ancTree");
//            File alignFile=new File(ptg.workDir.getParent()+File.separator+"6_leaves_test_set.aln");
//            ptg.percentPruning=1.0;
            //pplacer rRNA dataset
            
            //directory where all testes are done
//            File workDir=new File("/media/ben/STOCK/DATA/viromeplacer/accu_tests");
//            File treeFile=new File(workDir+File.separator+"RAxML_result.bv_refs_aln");
//            File alignFile=new File(workDir+File.separator+"bv_refs_aln_stripped_99.5.fasta");
//            //set directory associated to this ancTree
//            ptg.workDir=new File(workDir.getAbsolutePath()+File.separator+"pplacer_16s");


//            
//            //read ancTree
//            BufferedReader br=new BufferedReader(new FileReader(treeFile));
//            String line=null;
//            String treeString=null;
//            while((line=br.readLine())!=null) {treeString=line;}
//            PhyloTree ancTree = NewickReader.parseNewickTree2(treeString, true);
//            ancTree.initIndexes();
//            System.out.println("Original ancTree, # nodes: "+ancTree.getNodeCount());
//            //tree.displayTree();

                    
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
            //Ax: the pruned alignments
            //Tx: the pruned trees
            //Rx: the reads built from the pruned leaves
            //AND PREPARE THE AR direcoties FOR ALL minK/alpha COMBINATIONS (Dx)
            //i.e : build the extended trees
            //      build AR command, without execution
            //      prepare script which launch these using SGE (qsub)  
            //AND BUILD HMM directories and alignments for all pruning experiments
            //i.e. : build the hmm profile from the corresponding AxFile
            //       align all Rx reads to this profile
            //       convert the result from stockolm to fasta alignment
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
        File AxDir=new File(workDir+File.separator+"Ax");
        File TxDir=new File(workDir+File.separator+"Tx");
        File RxDir=new File(workDir+File.separator+"Rx");
        File GxDir=new File(workDir+File.separator+"Gx");
        AxDir.mkdir();
        TxDir.mkdir();
        RxDir.mkdir();
        GxDir.mkdir();
            
        //File where to save all model optimisations qsub commands
        File qSubOptimCommands=new File(workDir+File.separator+"Tx"+File.separator+"qsub_optim_commands");
        BufferedWriter bwOptim=new BufferedWriter(new FileWriter(qSubOptimCommands));
        
        //write in a binary file the pruned ancTree and the expected placement,
        //that is the branch b_new (and b_new_p if rerooting), see algo below.
        //expected placement is saved through an integer array of 1 or 2 elements
        //integer is the n of the node son of b_new
        File expectedPlacementsFile=new File(workDir.getAbsolutePath()+File.separator+"expected_placements.bin");
        ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(expectedPlacementsFile),4096));
        HashMap<Integer,Integer> NxIndex=new HashMap<>(); //map(120,it.e. nx120)=0; map(1142)=1 ; map(5454)=2 ...
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
        
        //all AR reconstructions and viromeplacer DBs will be in directory Dx
        File Dx=new File(workDir.getAbsolutePath()+File.separator+"Dx");
        Dx.mkdir();

        //HMM directory and qsub file
        //prepare files with all commands necessary for the HMM alignment build
        //1. commands for 1 alignment are saved as AxFile "set" of commands
        //in AxFile sh script.
        //2. all sh script calls are associated to qsub commands and put in
        //a single file
        StringBuilder qSubHMMCommands=new StringBuilder();
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
        for (int i = 0; i < prunedNodeIds.length; i++) {
            Integer nx_nodeId = prunedNodeIds[i];
            //System.out.println("--------------------------------------");
            //System.out.println("copying ancTree and alignment before pruning...");
            PhyloNode rootCopy=tree.getRoot().copy();
            PhyloTree treeCopy=new PhyloTree(new PhyloTreeModel(rootCopy),tree.isRooted(), false);
            //System.out.println("indexing ancTree ...");
            treeCopy.initIndexes();
            //some checkup about the original  copy
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
                Enumeration<PhyloNode> DFSenum = Nx.depthFirstEnumeration();
                //careful, ancTree topology change are reflected in the enumeration
                //so we need to remove node AFTER this this transversal postorder
                while (DFSenum.hasMoreElements()) {
                    nextNx=DFSenum.nextElement();
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
            
            

            
            File AxFile=new File(AxDir+File.separator+"A"+i+"_nx"+nx_nodeId+"_la"+Nx.getLabel()+".align");
            File TxFile=new File(TxDir+File.separator+"T"+i+"_nx"+nx_nodeId+"_la"+Nx.getLabel()+".tree");
            File GxFile=new File(GxDir+File.separator+"G"+i+"_nx"+nx_nodeId+"_la"+Nx.getLabel()+".fasta");
            
            // !!!!!!!!!!!!!!
            // This label is the identifier used to describe this particular pruning experiment !
            String experimentLabel=AxFile.getName().split("\\.align$")[0];
            // !!!!!!!!!!!!!!
            
            
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
            //save the pruned ancTree
            System.out.println("Indexing pruned tree");
            treeCopy.initIndexes(); //necessary to remove former nodes from the maps shortcuts
            System.out.println("pruned tree(treeCopy), #nodes :"+treeCopy.getNodeCount());
            System.out.println("pruned tree(treeCopy), #leaves:"+treeCopy.getLeavesCount());
            //it is necessary to keep the trees objects, if saved as
            //newick, then reloaded, the ancTree node ids will be different
            //and the expectedPlacement map will not match
            prunedTrees.add(treeCopy);
            System.out.println("Writing pruned tree newick");
            NewickWriter nw=new NewickWriter(TxFile);            
            nw.writeNewickTree(treeCopy, true, true, false, false);  //no internal node names if PAML !
            nw.close();
            
            //prepare bin data that will be used to compute trifurcations NDs
            //will not be computed if trifurcation  not asked for this Nx
            prunedTreesTrifurcations.add(new ArrayList<>());
            if (trifurcations && trifurcationsNxIndexes.contains(i)) {
                //generate trifurcations
                ArrayList<Integer> nodeIdsByDFS = treeCopy.getNodeIdsByDFS();
                int numberTrifurcationsExpected=treeCopy.getLeavesCount()-3;
                System.out.println("Will generate "+numberTrifurcationsExpected+" trifurcations. ");
                int counter=0;
                for (int it = 0; it < nodeIdsByDFS.size(); it++) {
                    Integer nodeId = nodeIdsByDFS.get(it);
                    if (!treeCopy.getById(nodeId).isLeaf()) {
                        //do not reroot on original root
                        if (treeCopy.getById(nodeId).isRoot()) {
                            System.out.println("Skipping node of original root.");
                            continue;
                        }
                        if ((counter%100)==0) {
                            System.out.println(counter+"/"+numberTrifurcationsExpected);
                        }
                        //copy tree
                        PhyloNode rootCopy2=treeCopy.getRoot().copy();
                        PhyloTree treeCopy2=new PhyloTree(new PhyloTreeModel(rootCopy2),treeCopy.isRooted(), false);
                        treeCopy2.initIndexes();
                        //reroot on nodeId
                        treeCopy2.rerootTree(treeCopy2.getById(nodeId), true);
                        //output as newick
                        File outname=new File(TxDir+File.separator+"trifu_"+counter+"_T"+i+"_nx"+nx_nodeId+"_la"+Nx.getLabel()+".tree");
                        BufferedWriter bwOut = Files.newBufferedWriter(outname.toPath());
                        new NewickWriter(bwOut).writeNewickTree(treeCopy2, true, true, false, false);
                        bwOut.close();
                        //save tree in experiment data
                        //it is important that this tree indexes are not reinitilaized
                        //to keep the correct same and be able to use the same Dtx
                        prunedTreesTrifurcations.get(prunedTreesTrifurcations.size()-1).add(treeCopy2);
                        
                        //discard copy
                        //treeCopy2=null; 
                        counter++;
                    }
                }
                System.out.println("# trifurcations generated: "+(counter+1));
            }
            
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
            
          
//            System.out.println("----------------");
//            System.out.println(nodeDistMatrix.toString());
//            System.out.println("----------------");
//            System.out.println(branchDistMatrix.toString());
//            System.out.println("----------------");

            ////////////////////////////////////////////////////////////////////
            //fourth, build Rx virtual read datasets from removedLeaves
                 
            //Save this pruned leafs in Gx, which can be used later for 
            //simulating even more reads
            FileWriter fwGx=new FileWriter(GxFile);
            for (Iterator<String> it = leavesRemoved.iterator(); it.hasNext();) {
                String next = it.next();
                fwGx.append(align.getFasta(next, false).getFormatedFasta());
                fwGx.append("\n");
            }
            fwGx.close();
                    
            String expString="R"+i+"_nx"+nx_nodeId+"_la"+Nx.getLabel();
            for (int j = 0; j < R.length; j++) {
                System.out.println("Preparing "+(j+1)+"th query read length in R="+Arrays.toString(R)+" ; Rsd="+Rsd);
                //prepare the corresponding output files
                String readFileString=expString+"_r"+R[j]+".fasta";
                File Rxj=new File(RxDir+File.separator+readFileString);
                //register it
                if (!readFiles.containsKey(experimentLabel))
                    readFiles.put(experimentLabel,new ArrayList<String>(R.length));
                readFiles.get(experimentLabel).add(Rxj.getName());
                       
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
                        int v_length = new Double(0.0+r+rand.nextGaussian()*Rsd).intValue();
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
            //fifth, build qsub commands to generate model parameters
            //for this particular pruning

            StringBuilder optimCommand=new StringBuilder();
            optimCommand.append(RAXMLExecutable.getAbsolutePath());
            optimCommand.append(" -f e");
            if (states instanceof AAStates) {
                optimCommand.append(" -m PROTGAMMAWAG -c 4 ");
            } else {
                optimCommand.append(" -m GTRGAMMA -c 4 ");
            };
            optimCommand.append(" -s ");
            optimCommand.append(AxFile.getAbsolutePath());
            optimCommand.append(" -t ");
            optimCommand.append(TxFile.getAbsolutePath());
            optimCommand.append(" -n OPTIM_");                    
            optimCommand.append(TxFile.getName());                    
            bwOptim.append("echo \""+optimCommand.toString()+"\" | qsub -N OPTIM_"+TxFile.getName()+" -wd "+TxDir.getAbsolutePath());
            bwOptim.newLine();
            
            
            
            ////////////////////////////////////////////////////////////////////
            //sixth, launch the generation of extended trees (used in AR) 
            //for this particular pruning
            
            
            //load data necessary to do the extended ancTree            
            String experimentAlignmentLabel=AxFile.getName().split("\\.align$")[0];
            
            System.out.println("PREPARING EXTENDED TREE/ALL DB_BUILDs FOR EXPERIMENT: "+experimentAlignmentLabel);

            //alignment filename is used as default workDir for this experiment
            File DxExpPath=new File(Dx+File.separator+experimentAlignmentLabel);
            System.out.println("Experiment work dir: "+DxExpPath.getAbsolutePath());
            if (workDir.canWrite())
                DxExpPath.mkdir();
            else {
                System.out.println("Cannot write in dir: "+DxExpPath.getAbsolutePath());
                System.exit(1);
            }
            
            //directory of extended ancTree
            File extendedTreeDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"extended_trees");
            extendedTreeDir.mkdir();
            
            File fileRelaxedAlignmentPhylip=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.phylip");
            File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
            File fileRelaxedAlignmentFasta=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.fasta");
            File fileRelaxedTreewithBL=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
//            File fileRelaxedTreeBinary=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree.bin");
            
            //build the extended trees first and save them in the extended_trees directory
            //these will be used in the AR below
            System.out.println("Create extended tree and save newick/fasta/phylip.");
            
            //do a copy or else, the extended ancTree will be saved in expected_placements.bin
            PhyloNode rootCopy2=treeCopy.getRoot().copy();
            PhyloTree treeCopy2=new PhyloTree(new PhyloTreeModel(rootCopy2),treeCopy.isRooted(), false);
            treeCopy2.initIndexes();
            System.out.println("treecopy2, # nodes: "+treeCopy2.getNodeCount());
            System.out.println("treecopy2, # leaves: "+treeCopy2.getLeavesCount());
            
            ExtendedTree extendedTree = buildExtendedTrees(extendedTreeDir, alignCopy, treeCopy2);         
            System.out.println("extendedTree, # nodes: "+extendedTree.getNodeCount());
            System.out.println("extendedTree, # leaves: "+extendedTree.getLeavesCount());       
            
            ////////////////////////////////////////////////////////////////////
            //seventh, make subdirectory corresponding to  minK/alpha combinations
            //this will be the directory in which viromeplacer will work later
            
            //prepare the directory for AR based on extended ancTree 
            //create correpsonding directory        
            File ARDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"AR");
            ARDir.mkdir();
            
            for (int current_k = minK; current_k < maxK+kIncrement; current_k+=kIncrement) {
                for (float current_alpha = minFactor; current_alpha < maxFactor+factorIncrement; current_alpha+=factorIncrement) {

                    //System.out.print("; minK="+current_k+" a="+nf.format(current_alpha));

                    //make subdirectory corresponding to this minK/alpha combination
                    File combDir=new File(DxExpPath.getAbsolutePath()+File.separator+"k"+current_k+"_a"+nf.format(current_alpha));
                    boolean dirOK=combDir.mkdir();
                    
                    //do AxFile subdir AR/ in each minK/alpha directory
                    File combDirAR=new File(combDir.getAbsolutePath()+File.separator+"AR");
                    combDirAR.mkdir();
                    File combDirExtendedTree=new File(combDir.getAbsolutePath()+File.separator+"extended_trees");
                    combDirExtendedTree.mkdir();
                    
                    //make symbolic link to rst file which is in experiment directory (same AR for all minK/alpha combinations)
                    File originalRST=new File(ARDir.getAbsolutePath()+File.separator+"rst");
                    originalRST.createNewFile();
                    File linkToRST=new File(combDirAR.getAbsolutePath()+File.separator+"rst");
                    if (!linkToRST.exists())
                        Files.createSymbolicLink(linkToRST.toPath(), originalRST.toPath());
                    
                    //make symbolic link to phyML AR outputs files which are in experiment directory (same AR for all minK/alpha combinations)
                    File originalAncestralSeq=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_seq.txt");
                    originalAncestralSeq.createNewFile();
                    File linkToAncestralSeq=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_seq.txt");
                    if (!linkToAncestralSeq.exists())
                        Files.createSymbolicLink(linkToAncestralSeq.toPath(), originalAncestralSeq.toPath());
                    
                    File originalAncestralTree=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_tree.txt");
                    originalAncestralTree.createNewFile();
                    File linkToAncestralTree=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_tree.txt");
                    if (!linkToAncestralTree.exists())
                        Files.createSymbolicLink(linkToAncestralTree.toPath(), originalAncestralTree.toPath());
                    
                    File originalPhymlStats=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_stats.txt");
                    originalPhymlStats.createNewFile();
                    File linkToOriginalPhymlStats=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_stats.txt");
                    if (!linkToOriginalPhymlStats.exists())
                        Files.createSymbolicLink(linkToOriginalPhymlStats.toPath(), originalPhymlStats.toPath());
                    
                    File originalPhymlTree=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_tree.txt");
                    originalPhymlTree.createNewFile();
                    File linkToOriginalPhymlTree=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_tree.txt");
                    if (!linkToOriginalPhymlTree.exists())
                        Files.createSymbolicLink(linkToOriginalPhymlTree.toPath(), originalPhymlTree.toPath());
                    
                    //do AxFile subdir extended_trees/ in each minK/alpha directory
                    //make symbolic link to correpsonding trees/alignments
                    //which are in experiment directory (same trees/aligns for all minK/alpha combinations)
                    //this is necessary because these trees are required by db_build
                    File linkToTree=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
                    File linkToTreeNoInter=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
                    File linkToAlignmentFasta=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_align.fasta");
                    File linkToAlignmentPhylip=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_align.phylip");
//                    File linkToRelaxedTreeBinary=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree.bin");
                    
                    if (!linkToTree.exists())
                        Files.createSymbolicLink(linkToTree.toPath(), fileRelaxedTreewithBL.toPath());       
                    if (!linkToTreeNoInter.exists())
                        Files.createSymbolicLink(linkToTreeNoInter.toPath(), fileRelaxedTreewithBLNoInternalNodeLabels.toPath());  
                    if (!linkToAlignmentFasta.exists())
                        Files.createSymbolicLink(linkToAlignmentFasta.toPath(), fileRelaxedAlignmentFasta.toPath()); 
                    if (!linkToAlignmentPhylip.exists())
                        Files.createSymbolicLink(linkToAlignmentPhylip.toPath(), fileRelaxedAlignmentPhylip.toPath());
//                    if (!linkToRelaxedTreeBinary.exists()) 
//                        Files.createSymbolicLink(linkToRelaxedTreeBinary.toPath(), fileRelaxedTreeBinary.toPath());
                    
                }
            }
            
            
            ////////////////////////////////////////////////////////////////////
            //seventh-BIS (optionnal) trifurcations versions
            //make subdirectory corresponding to  minK/alpha combinations
            //this will be the directory in which viromeplacer will work later

            //if trifurcations, complete with as many xtended_trees variants
            //directory naming will be extended_trees_trifuX
            if (trifurcations && trifurcationsNxIndexes.contains(i)) {
                //generate trifurcations
                ArrayList<Integer> nodeIdsByDFS = treeCopy.getNodeIdsByDFS();
                int numberTrifurcationsExpected=treeCopy.getLeavesCount()-3;
                System.out.println("creating directories and links extended_trees_trifuX & AR_trifuX");
                int counter=0;
                for (int it = 0; it < nodeIdsByDFS.size(); it++) {
                    Integer nodeId = nodeIdsByDFS.get(it);
                    if (!treeCopy.getById(nodeId).isLeaf()) {
                        //do not reroot on original root
                        if (treeCopy.getById(nodeId).isRoot()) {
                            System.out.println("Skipping node of original root.");
                            continue;
                        }
                        if ((counter%100)==0) {
                            System.out.println(counter+"/"+numberTrifurcationsExpected);
                        }
                        //if (counter>20) {break;}
                        
                        //directory of AR 
                        ARDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"AR_trifu"+counter);
                        ARDir.mkdir();
                        
                        //directory of extended ancTree
                        extendedTreeDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"extended_trees_trifu"+counter);
                        extendedTreeDir.mkdir();

                        fileRelaxedAlignmentPhylip=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.phylip");
                        fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
                        fileRelaxedAlignmentFasta=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.fasta");
                        fileRelaxedTreewithBL=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
        //              File fileRelaxedTreeBinary=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree.bin");

                        //build the extended trees first and save them in the extended_trees directory
                        //these will be used in the AR below
                        //do a copy or else, the extended ancTree will be saved in expected_placements.bin
                        rootCopy2=prunedTreesTrifurcations.get(prunedTreesTrifurcations.size()-1).get(counter).getRoot().copy();
                        treeCopy2=new PhyloTree(new PhyloTreeModel(rootCopy2),prunedTreesTrifurcations.get(prunedTreesTrifurcations.size()-1).get(counter).isRooted(), false);
                        treeCopy2.initIndexes();
                        
                        //load correct pruned alignment
                        FASTAPointer fp=new FASTAPointer(AxFile, false);
                        Fasta fasta=null;
                        ArrayList<Fasta> fastas=new ArrayList<>();
                        while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                            fastas.add(fasta);
                        }
                        Alignment alignCopy2=new Alignment(s,fastas);
                        fp.closePointer();
                        
                        
                        extendedTree = buildExtendedTrees(extendedTreeDir, alignCopy2, treeCopy2);
                        
                        extendedTree=null;
                        treeCopy2=null;
                        alignCopy2=null;
                        
                        
                        for (int current_k = minK; current_k < maxK+kIncrement; current_k+=kIncrement) {
                            for (float current_alpha = minFactor; current_alpha < maxFactor+factorIncrement; current_alpha+=factorIncrement) {

                                //make subdirectory corresponding to this minK/alpha combination
                                File combDir=new File(DxExpPath.getAbsolutePath()+File.separator+"k"+current_k+"_a"+nf.format(current_alpha));

                                //do AxFile subdir AR/ in each minK/alpha directory
                                File combDirAR=new File(combDir.getAbsolutePath()+File.separator+"AR_trifu"+counter);
                                combDirAR.mkdir();
                                File combDirExtendedTree=new File(combDir.getAbsolutePath()+File.separator+"extended_trees_trifu"+counter);
                                combDirExtendedTree.mkdir();

                                //make symbolic link to rst file which is in experiment directory (same AR for all minK/alpha combinations)
                                File originalRST=new File(ARDir.getAbsolutePath()+File.separator+"rst");
                                originalRST.createNewFile();
                                File linkToRST=new File(combDirAR.getAbsolutePath()+File.separator+"rst");
                                if (!linkToRST.exists())
                                    Files.createSymbolicLink(linkToRST.toPath(), originalRST.toPath());

                                //make symbolic link to phyML AR outputs files which are in experiment directory (same AR for all minK/alpha combinations)
                                File originalAncestralSeq=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_seq.txt");
                                originalAncestralSeq.createNewFile();
                                File linkToAncestralSeq=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_seq.txt");
                                if (!linkToAncestralSeq.exists())
                                    Files.createSymbolicLink(linkToAncestralSeq.toPath(), originalAncestralSeq.toPath());

                                File originalAncestralTree=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_tree.txt");
                                originalAncestralTree.createNewFile();
                                File linkToAncestralTree=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_ancestral_tree.txt");
                                if (!linkToAncestralTree.exists())
                                    Files.createSymbolicLink(linkToAncestralTree.toPath(), originalAncestralTree.toPath());

                                File originalPhymlStats=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_stats.txt");
                                originalPhymlStats.createNewFile();
                                File linkToOriginalPhymlStats=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_stats.txt");
                                if (!linkToOriginalPhymlStats.exists())
                                    Files.createSymbolicLink(linkToOriginalPhymlStats.toPath(), originalPhymlStats.toPath());

                                File originalPhymlTree=new File(ARDir.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_tree.txt");
                                originalPhymlTree.createNewFile();
                                File linkToOriginalPhymlTree=new File(combDirAR.getAbsolutePath()+File.separator+"extended_align.phylip_phyml_tree.txt");
                                if (!linkToOriginalPhymlTree.exists())
                                    Files.createSymbolicLink(linkToOriginalPhymlTree.toPath(), originalPhymlTree.toPath());

                                //do AxFile subdir extended_trees/ in each minK/alpha directory
                                //make symbolic link to correpsonding trees/alignments
                                //which are in experiment directory (same trees/aligns for all minK/alpha combinations)
                                //this is necessary because these trees are required by db_build
                                File linkToTree=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
                                File linkToTreeNoInter=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
                                File linkToAlignmentFasta=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_align.fasta");
                                File linkToAlignmentPhylip=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_align.phylip");
            //                    File linkToRelaxedTreeBinary=new File(combDirExtendedTree.getAbsolutePath()+File.separator+"extended_tree.bin");

                                if (!linkToTree.exists())
                                    Files.createSymbolicLink(linkToTree.toPath(), fileRelaxedTreewithBL.toPath());       
                                if (!linkToTreeNoInter.exists())
                                    Files.createSymbolicLink(linkToTreeNoInter.toPath(), fileRelaxedTreewithBLNoInternalNodeLabels.toPath());  
                                if (!linkToAlignmentFasta.exists())
                                    Files.createSymbolicLink(linkToAlignmentFasta.toPath(), fileRelaxedAlignmentFasta.toPath()); 
                                if (!linkToAlignmentPhylip.exists())
                                    Files.createSymbolicLink(linkToAlignmentPhylip.toPath(), fileRelaxedAlignmentPhylip.toPath());
            //                    if (!linkToRelaxedTreeBinary.exists()) 
            //                        Files.createSymbolicLink(linkToRelaxedTreeBinary.toPath(), fileRelaxedTreeBinary.toPath());
      
                            }//end for alpha
                        }//end for k

                        counter++;    
                        
                    }//end if isLeaf
                }//end for DFSids 
            }//end if trifurcation

                           
            ////////////////////////////////////////////////////////////////////
            //seventh, build HMM command matching this pruning experiment
            //
            
            //alignment filename is used as default workDir for this experiment
            File HMMxAxDir=new File(HMMxDir.getAbsolutePath()+File.separator+experimentAlignmentLabel);
            System.out.println("HMM Experiment work dir: "+HMMxAxDir.getAbsolutePath());
            if (workDir.canWrite())
                HMMxAxDir.mkdir();
            else {
                System.out.println("Cannot write in dir: "+HMMxAxDir.getAbsolutePath());
                System.exit(1);
            }
            
            //all files used in this HMM alignment
            File hmmOuput=new File(HMMxAxDir.getAbsolutePath()+File.separator+experimentAlignmentLabel+".hmm");
            
            //build the hmmbuild command. 1 hmm per experiment
            StringBuilder commandSet=new StringBuilder();
            //commandSet.append("#!/bin/sh\n");
            commandSet.append(HMMBUILDPath.getAbsolutePath()+" ");
            if (proteinAnalysis) {
                commandSet.append("--amino ");
            } else {
                commandSet.append("--dna ");
            }
            commandSet.append(hmmOuput.getAbsolutePath()+" "+AxFile.getAbsolutePath()+"\n");
            
            //build alignments of reads
            List<String> reads=readFiles.get(experimentAlignmentLabel);
            //System.out.println("Rx files: "+readFiles.get(experimentLabel));
            for (int j=0;j<reads.size();j++) {
                File alnOutput=new File(HMMxAxDir.getAbsolutePath()+File.separator+reads.get(j).split("\\.fasta$")[0]+".aln.psiblast");
                commandSet.append(HMMALIGNPath.getAbsolutePath()+" --outformat PSIBLAST ");
                if (proteinAnalysis) {
                    commandSet.append("--amino ");
                } else {
                    commandSet.append("--dna ");
                }
                commandSet.append("--mapali "+AxFile.getAbsolutePath()+" -o "+alnOutput.getAbsolutePath()
                                    + " "+hmmOuput+" "+RxDir.getAbsolutePath()+File.separator+reads.get(j)
                                    + "\n");     
                //conversion to fasta
                commandSet.append(  pyScript.getAbsolutePath()+" "+
                                    alnOutput.getAbsolutePath()+" "+
                                    alnOutput.getParentFile().getAbsolutePath()+File.separator+alnOutput.getName().split("\\.aln.psiblast$")[0]+".aln.fasta" +
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
            qSubHMMCommands.append("echo \""+shScript.getAbsolutePath()+"\" |  qsub -N HMM_"+experimentAlignmentLabel+
                                " -wd "+HMMxAxDir.getAbsolutePath()+"\n");
            
            
            
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
        //give execution rights to the commands file
        Files.setPosixFilePermissions(qSubOptimCommands.toPath(), perms);
        bwOptim.close();
        
//        for (int it = 0; it < prunedTrees.size(); it++) {
//            PhyloTree get = prunedTrees.get(it);
//            System.out.println(it+"th ancTree size test:"+get.getNodeCount());
//        }
        
        //after all placements, save expected placements in binary file
        oos.writeObject(NxIndex);
        oos.writeObject(expectedPlacements);
        oos.writeObject(prunedTrees);
        oos.writeObject(prunedTreesTrifurcations);
            
        //save qsub HMM commands
        File qsubScript=new File(HMMxDir.getAbsolutePath()+File.separator+"qsub_hmm_commands");
        BufferedWriter bwQsubScript=new BufferedWriter(new FileWriter(qsubScript));
        bwQsubScript.append(qSubHMMCommands.toString());
        bwQsubScript.close();
        Files.setPosixFilePermissions(qsubScript.toPath(), perms);
        
        
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
    
 
    /**
     * build extendedTree from given PhyloTree and Alignment, 
     * save it as newick, phylip, fasta and binary files in extendedTreePath
     * @param extendedTreePath
     * @param align
     * @param tree 
     */
    private ExtendedTree buildExtendedTrees(File extendedTreePath,Alignment align, PhyloTree tree) {
        
        ExtendedTree extendedTree = null;

        File fileRelaxedAlignmentFasta=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_align.fasta");
        File fileRelaxedAlignmentPhylip=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_align.phylip");
        File fileRelaxedTreewithBL=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
        File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");;
        //File fileRelaxedTreeBinary=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_tree.bin");
        //File idsMappings=new File(extendedTreePath.getAbsolutePath()+File.separator+"extended_tree_node_mapping_from_prunedTreeGenerator.tsv");
        
        try {
            //note, we read again the ancTree to build AxFile new PhyloTree object
            //this is necessary as its TreeModel is directly modified
            //at instanciation of ExtendedTree
            extendedTree=new ExtendedTree(tree,minBranchLength,branchPerEdge);                    
            extendedTree.initIndexes(); // don'TxFile forget to reinit indexes !!!
            ArrayList<PhyloNode> listOfNewFakeLeaves = extendedTree.getFakeLeaves();
            //System.out.println("RelaxedTree contains "+extendedTree.getLeavesCount()+ " leaves");
            //System.out.println("RelaxedTree contains "+extendedTree.getFakeLeaves().size()+ " FAKE_X new leaves");
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
            //System.out.println("Write extended alignment (fasta): "+fileRelaxedAlignmentFasta.getAbsolutePath());
            align.writeAlignmentAsFasta(fileRelaxedAlignmentFasta);
            //System.out.println("Write extended alignment (phylip): "+fileRelaxedAlignmentPhylip.getAbsolutePath());
            align.writeAlignmentAsPhylip(fileRelaxedAlignmentPhylip);
            align=null;
            //write extended trees
            //System.out.println("Write extended newick ancTree: "+fileRelaxedTreewithBL.getAbsolutePath());
            NewickWriter nw=new NewickWriter(fileRelaxedTreewithBL);
            nw.writeNewickTree(extendedTree, true, true, false, false);
            nw.close();
            //write version without internal nodes labels
            //System.out.println("Write extended newick ancTree with branch length: "+fileRelaxedTreewithBLNoInternalNodeLabels.getAbsolutePath());
            nw=new NewickWriter(fileRelaxedTreewithBLNoInternalNodeLabels);
            nw.writeNewickTree(extendedTree, true, false, false, false);
            nw.close();
            //save this extendedTree as AxFile binary
//            FileOutputStream fos = new FileOutputStream(fileRelaxedTreeBinary);
//            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
//            System.out.println("Storing binary version of ExtendedTree: "+fileRelaxedTreeBinary.getAbsolutePath());
//            oos.writeObject(extendedTree);
//            oos.close();
            //finally, for debugging, output the ids mappings
//            FileWriter fwRxj=new FileWriter(idsMappings);
//            fwRxj.append("original_id\toriginal_name\textended_id\textended_name");
//            LinkedHashMap<Integer, Integer> fakeNodeMapping = extendedTree.getFakeNodeMapping();
//            for (Iterator<Integer> iterator = fakeNodeMapping.keySet().iterator(); iterator.hasNext();) {
//                Integer next = iterator.next();
//                fwRxj.append("\n");
//                fwRxj.append(fakeNodeMapping.get(next)+"\t"+ancTree.getById(fakeNodeMapping.get(next)).getLabel()+"\t");
//                fwRxj.append(next+"\t"+extendedTree.getById(next).getLabel());
//            }
//            fwRxj.close();  
        } catch (IOException ex) {
            ex.printStackTrace();
            System.out.println("Error raised from extended tree reconstruction!");
        }
        return extendedTree;
    }
    
}
