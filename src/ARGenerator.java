
import alignement.Alignment;
import core.AAStates;
import core.DNAStatesShifted;
import core.States;
import inputs.ARProcessLauncher;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.file.Files;
import java.nio.file.attribute.PosixFilePermission;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import main_v2.ArgumentsParser_v2;
import models.EvolModel;
import tree.NewickReader;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.PhyloTreeModel;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class ARGenerator {
//pseudorandom generator
    Long seed=new Long(1);
    Random rand = new Random(seed);
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");
    
    //default dataset
    File alignFile=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set.aln");
    File treeFile=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set.tree");
    
    //nodes selected for pruning at this launch
    Integer[] prunedNodeIds=null;
    //nodes effectively pruned because didn'TxFile produce tre of less than 4 leaves
    ArrayList<Integer> actuallyPrunedNodeIdsIndexes=new ArrayList<>();
    
    //list of new files
    public HashMap<String,List<String>> readFiles=new HashMap<>(); //map of reads, map(A0_nx4_la)=[A0_nx4_la_r150,A0_nx4_la_r300]
    public ArrayList<PhyloTree> prunedTrees=new ArrayList<>(); //list of pruned Trees
    
    //pruning fraction
    //double percentPruning=0.5; //10%
    int pruningCount=100;

    
    //set if analysis is protein or DNA/RNA
    boolean proteinAnalysis=false;
    
    //Generate all possible trifurcations for each pruning
    //build as many databases and placement commands.
    //The same parameters and alignment file is used.
    //trifurcations follow naming convention :
    //Tx_trifuXX_nxXX_laXX
    boolean trifurcations=false;
    ArrayList<Integer> trifurcationsNxIndexes=null;

    
    //AR program
    File ARExecutablePath=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/paml4.9b_hacked/bin/baseml");

    
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
    
    //load the expected placement
    HashMap<Integer,Integer> NxIndex =null;
    ArrayList<PhyloTree> experimentTrees = null;
    ArrayList<ArrayList<PhyloTree>> experimentTreesTrifurcations =null;
    
        
    public static void main(String[] args) {
        
        System.out.println("ARGS: workDir ARBinaries align tree percentPruning(float) [nucl=0|prot=1] [trifurcations:-1=no|1,45,48=list of Nx to test] ");
        
        System.out.println("Command: "+Arrays.toString(args).replaceAll(",", " "));

        
        try {
            //launch
            ARGenerator arg=new ARGenerator();
            
            //PREPARE DATA SOURCES
            ///////////////////////////////////////////////////
            if (args.length>0) {
                arg.workDir=new File(args[0]);
                arg.ARExecutablePath=new File(args[1]);
                arg.alignFile=new File(args[2]);
                arg.treeFile=new File(args[3]);
                //ptg.percentPruning=Double.parseDouble(args[5]);
                arg.pruningCount=Integer.parseInt(args[4]);
                if (args.length>5) {
                    int protein=Integer.parseInt(args[5]);
                    arg.proteinAnalysis=(protein>0);
                if (!args[6].equals("-1")) {
                        arg.trifurcations=true;
                        String[] trifuIndexes=args[6].split(",");
                        arg.trifurcationsNxIndexes=new ArrayList<>(trifuIndexes.length);
                        for (String trifuIndexe : trifuIndexes) {
                            arg.trifurcationsNxIndexes.add(Integer.valueOf(trifuIndexe));
                        }
                    }
                }
            }
            
            System.out.println("workDir: "+arg.workDir);
            System.out.println("ARExecutable: "+arg.ARExecutablePath);
            System.out.println("alignFile: "+arg.alignFile);
            System.out.println("treeFile: "+arg.treeFile);
            //System.out.println("percentPruning: "+arg.percentPruning);
            System.out.println("pruningCount: "+arg.pruningCount);
            System.out.println("protein:"+arg.proteinAnalysis);
            System.out.println("trifurcations:"+arg.trifurcations);
            System.out.println("trifurcations Nx tested: "+arg.trifurcationsNxIndexes);

            //TEST ZONE
            
            //very small default dataset for debugging
//            File treeFile=new File(arg.workDir.getParent()+File.separator+"6_leaves_test_set.ancTree");
//            File alignFile=new File(arg.workDir.getParent()+File.separator+"6_leaves_test_set.aln");
//            arg.percentPruning=1.0;
            //pplacer rRNA dataset
            
            //directory where all testes are done
//            File workDir=new File("/media/ben/STOCK/DATA/viromeplacer/accu_tests");
//            File treeFile=new File(workDir+File.separator+"RAxML_result.bv_refs_aln");
//            File alignFile=new File(workDir+File.separator+"bv_refs_aln_stripped_99.5.fasta");
//            //set directory associated to this ancTree
//            arg.workDir=new File(workDir.getAbsolutePath()+File.separator+"pplacer_16s");


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
            States s=new DNAStatesShifted();
            if (arg.proteinAnalysis) {
                s=new AAStates(true);
            }
            FASTAPointer fp=new FASTAPointer(arg.alignFile, false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
            }
            Alignment align=new Alignment(s,fastas);
            System.out.println(align.describeAlignment(false));
            fp.closePointer();
            //load ancTree
            BufferedReader br=new BufferedReader(new FileReader(arg.treeFile));
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
            
            //if trifurcations are tested, load pruned trees
            if (arg.trifurcations) {
                File expPLaceFile=new File(arg.workDir.getAbsolutePath()+File.separator+"expected_placements.bin");
                try {
                    ObjectInputStream ois=null;
                    ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(expPLaceFile),4096));
                    System.out.println("Loading NxIndex");
                    arg.NxIndex = (HashMap<Integer,Integer>)ois.readObject();
                    System.out.println("Loading expected placements");
                    ArrayList<ArrayList<Integer>> expectedPlacementsNodeIds = (ArrayList<ArrayList<Integer>>)ois.readObject();
                    expectedPlacementsNodeIds=null;
                    System.out.println("Loading trees");
                    arg.experimentTrees = (ArrayList<PhyloTree>)ois.readObject();
                    System.out.println("Loading trees trifurcations");
                    arg.experimentTreesTrifurcations = (ArrayList<ArrayList<PhyloTree>>)ois.readObject(); 
                } catch (ClassNotFoundException ex) {
                    Logger.getLogger(ARGenerator.class.getName()).log(Level.SEVERE, null, ex);
                }
                if (arg.experimentTreesTrifurcations==null) {
                    System.out.println("Could not load trifurcation collection from expected_placements.bin");
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
            if (!arg.workDir.exists()) {
                arg.workDir.mkdir();
            }
            arg.generatePrunedTrees(arg.workDir,tree,align);
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

        
        //prepare the output directories
        File AxDir=new File(workDir+File.separator+"Ax");
        File TxDir=new File(workDir+File.separator+"Tx");

        //all AR reconstructions and viromeplacer DBs will be in directory Dx
        File Dx=new File(workDir.getAbsolutePath()+File.separator+"Dx");
        //File where to save all AR qsub commands
        File qSubARCommands=new File(workDir+File.separator+"Dx"+File.separator+"qsub_AR_commands");
        BufferedWriter bwAR=new BufferedWriter(new FileWriter(qSubARCommands));
        File qSubARCommandsTrifu=new File(workDir+File.separator+"Dx"+File.separator+"qsub_AR_commands_trifu");
        BufferedWriter bwARTrifu=new BufferedWriter(new FileWriter(qSubARCommandsTrifu));
        
        ////////////////////////////////////////////////////////////////////
        //First, reload correct prunings
        
        for (int i = 0; i < prunedNodeIds.length; i++) {
            Integer nx_nodeId = prunedNodeIds[i];            
            //current root Nx defining the pruned clade
            PhyloNode Nx= tree.getById(nx_nodeId);
            System.out.println("--------------------------------------");
            System.out.println("selected Nx: x="+i+"  node:"+Nx);
            if (Nx.isRoot()) {
                System.out.println("SKIPPED: this node is root (id="+nx_nodeId+"), no pruning.");
                continue;
            }
                        
            File AxFile=new File(AxDir+File.separator+"A"+i+"_nx"+nx_nodeId+"_la"+Nx.getLabel()+".align");
            File TxFile=new File(TxDir+File.separator+"T"+i+"_nx"+nx_nodeId+"_la"+Nx.getLabel()+".tree");

            ////////////////////////////////////////////////////////////////////
            //Second, launch the generation AR commands
            //for this particular pruning
            
            //load data necessary to do the extended ancTree            
            String experimentAlignmentLabel=AxFile.getName().split("\\.align$")[0];
            
            System.out.println("PREPARING AR FOR EXPERIMENT: "+experimentAlignmentLabel);

            //alignment filename is used as default workDir for this experiment
            File DxExpPath=new File(Dx+File.separator+experimentAlignmentLabel);
            
            //check if directory exists, if not it was a case with less than 3 leaves
            if (!DxExpPath.exists()) {
                System.out.println("SKIPPED: This directory do not exists (pruning results to a 3 with less than 3 leaves !)");              
                continue;
            }
            
            System.out.println("Experiment work dir: "+DxExpPath.getAbsolutePath());

            
            //load directory and files of extended ancTree (created at previous step)
            File extendedTreeDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"extended_trees");            
            File fileRelaxedAlignmentPhylip=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.phylip");
            File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
            File fileRelaxedAlignmentFasta=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.fasta");
            File fileRelaxedTreewithBL=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
//            File fileRelaxedTreeBinary=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree.bin");  

            //load the directory for AR based on extended ancTree  (created at previous step)
            File ARDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"AR");
            
            ////////////////////////////////////////////////////////////////////
            //third, build the AR commands themselves
            
            //load value of alpha from optimisation files
            File optimFile=new File(TxFile.getParentFile().getAbsolutePath()+File.separator+"RAxML_info.OPTIM_"+TxFile.getName());
            BufferedReader brOptimFile = Files.newBufferedReader(optimFile.toPath());
            String line=null;
            float a=1.0f;
            while ((line=brOptimFile.readLine())!=null) {      
                if (line.startsWith("alpha: ")) {
                    a=Float.parseFloat(line.split(" ")[1]);
                    
                    break;
                }
            }
            brOptimFile.close();
            System.out.println("Found in optim file: alpha="+a);
            
            System.out.println("Build AR commands in: "+ARDir.getAbsolutePath());
            ARProcessLauncher arpl=null;
            StringBuilder ARCommand = new StringBuilder("");
            if (proteinAnalysis) {
                EvolModel model=new EvolModel(ArgumentsParser_v2.STATES_PROTEIN, "LG", a, 4);
                arpl=new ARProcessLauncher(ARExecutablePath,false,new AAStates(false),model,null);
            } else {
                EvolModel model=new EvolModel(ArgumentsParser_v2.STATES_DNA, "GTR", a, 4);
                arpl=new ARProcessLauncher(ARExecutablePath,false,new DNAStatesShifted(),model,null);
            }
            ARCommand.append( arpl.prepareAR(ARDir, fileRelaxedAlignmentPhylip, fileRelaxedTreewithBLNoInternalNodeLabels) );
            
            if (ARExecutablePath.getName().contains("phyml")) {
                //add move of the 4 files to the AR directory, as phyml write 
                //outputs near its input alignment file
                //phyml is written all data files near the input aignment file...
                //we move them to the AR directory
                //files are:
                // 1. alignName_phyml_ancestral_seq.txt         (used)
                // 2. alignName_phyml_stats.txt                 (unused)
                // 3. alignName_phyml_ancestral_tree.txt        (used)
                // 4. alignName_phyml_tree.txt                  (unused)
                File stats=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_stats.txt");
                File ancTree=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_ancestral_tree.txt");
                File seq=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_ancestral_seq.txt");
                File oriTree=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_tree.txt");
                //rename to these files
                File statsNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_stats.txt");
                File ancTreeNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_ancestral_tree.txt");
                File seqNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_ancestral_seq.txt");
                File oriTreeNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_tree.txt");
                //add the mv commands
                ARCommand.append(" ; mv "+stats.getAbsolutePath()+" "+statsNew.getAbsolutePath());
                ARCommand.append(" ; mv "+ancTree.getAbsolutePath()+" "+ancTreeNew.getAbsolutePath());
                ARCommand.append(" ; mv "+seq.getAbsolutePath()+" "+seqNew.getAbsolutePath());
                ARCommand.append(" ; mv "+oriTree.getAbsolutePath()+" "+oriTreeNew.getAbsolutePath());
 
            } 
            
            //prepare corresponding AR commands (with PAML)
            bwAR.append("echo \""+ARCommand+"\" | qsub -N AR_"+DxExpPath.getName()+" -wd "+ARDir.getAbsolutePath());
            bwAR.newLine();
            
            
            ////////////////////////////////////////////////////////////////////
            //third-BIS (optionnal), build the AR commands for trifurcations
            //File where to save all AR qsub commands for trifurcations

            ARCommand = new StringBuilder("");
            
            //if trifurcations, complete with as many AR variants
            if (trifurcations && trifurcationsNxIndexes.contains(i)) {
                System.out.println("creating AR commands for trifurcations using extended_trees_trifuX & AR_trifuX");
                
                //load correct pruned tree
                PhyloTree prunedTree=experimentTrees.get(NxIndex.get(nx_nodeId));
                System.out.println("pruned tree, #nodes :"+prunedTree.getNodeCount());
                System.out.println("pruned tree, #leaves:"+prunedTree.getLeavesCount());
                
                //generate trifurcations
                ArrayList<Integer> nodeIdsByDFS = prunedTree.getNodeIdsByDFS();
                int numberTrifurcationsExpected=prunedTree.getLeavesCount()-3;
                int counter=0;
                for (int it = 0; it < nodeIdsByDFS.size(); it++) {
                    Integer nodeId = nodeIdsByDFS.get(it);
                    if (!prunedTree.getById(nodeId).isLeaf()) {
                        
                        //do not reroot on original root
                        if (prunedTree.getById(nodeId).isRoot()) {
                            System.out.println("Skipping node of original root.");
                            continue;
                        }
                        if ((counter%100)==0) {
                            System.out.println(counter+"/"+numberTrifurcationsExpected);
                        }
                        
                        //load directory and files of extended ancTree (created at previous step)
                        extendedTreeDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"extended_trees_trifu"+counter);            
                        fileRelaxedAlignmentPhylip=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.phylip");
                        fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
                        fileRelaxedAlignmentFasta=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_align.fasta");
                        fileRelaxedTreewithBL=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
            //            File fileRelaxedTreeBinary=new File(extendedTreeDir.getAbsolutePath()+File.separator+"extended_tree.bin");  

                        //load the directory for AR based on extended ancTree  (created at previous step)
                        ARDir=new File(workDir.getAbsolutePath()+File.separator+"Dx"+File.separator+experimentAlignmentLabel+File.separator+"AR_trifu"+counter);
                        ARCommand = new StringBuilder("");
                        ARCommand.append( arpl.prepareAR(ARDir, fileRelaxedAlignmentPhylip, fileRelaxedTreewithBLNoInternalNodeLabels) );

                        if (ARExecutablePath.getName().contains("phyml")) {
                            //add move of the 4 files to the AR directory, as phyml write 
                            //outputs near its input alignment file
                            //phyml is written all data files near the input aignment file...
                            //we move them to the AR directory
                            //files are:
                            // 1. alignName_phyml_ancestral_seq.txt         (used)
                            // 2. alignName_phyml_stats.txt                 (unused)
                            // 3. alignName_phyml_ancestral_tree.txt        (used)
                            // 4. alignName_phyml_tree.txt                  (unused)
                            File stats=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_stats.txt");
                            File ancTree=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_ancestral_tree.txt");
                            File seq=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_ancestral_seq.txt");
                            File oriTree=new File(fileRelaxedAlignmentPhylip.getAbsolutePath()+"_phyml_tree.txt");
                            //rename to these files
                            File statsNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_stats.txt");
                            File ancTreeNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_ancestral_tree.txt");
                            File seqNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_ancestral_seq.txt");
                            File oriTreeNew=new File(fileRelaxedAlignmentPhylip.getParent().replace("/extended_trees", "/AR")+File.separator+fileRelaxedAlignmentPhylip.getName()+"_phyml_tree.txt");
                            //add the mv commands
                            ARCommand.append(" ; mv "+stats.getAbsolutePath()+" "+statsNew.getAbsolutePath());
                            ARCommand.append(" ; mv "+ancTree.getAbsolutePath()+" "+ancTreeNew.getAbsolutePath());
                            ARCommand.append(" ; mv "+seq.getAbsolutePath()+" "+seqNew.getAbsolutePath());
                            ARCommand.append(" ; mv "+oriTree.getAbsolutePath()+" "+oriTreeNew.getAbsolutePath());

                        } 

                        //prepare corresponding AR commands (with PAML)
                        bwARTrifu.append("echo \""+ARCommand+"\" | qsub -N AR_"+DxExpPath.getName()+"_trifu"+counter+" -wd "+ARDir.getAbsolutePath());
                        bwARTrifu.newLine();
                        
                        
                        
                        counter++;
                    }
                }//enf for DFS nodeIds
            }//enf if trifurcation
            
            
            ////////////////////////////////////////////////////////////////////
            //closes everything
            
            // a garbage call to clean a bit before next pruning
            System.gc();
            
        } //END of nx_nodeId loop
        

        //give execution rights to the commands file
        Files.setPosixFilePermissions(qSubARCommands.toPath(), perms);
        bwAR.close();
        //give execution rights to the commands file
        Files.setPosixFilePermissions(qSubARCommandsTrifu.toPath(), perms);
        bwARTrifu.close();

         
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