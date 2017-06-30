
import alignement.Alignment;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import main_v2.ARProcessLauncher;
import main_v2.Main_v2;
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
 *
 * @author ben
 */
public class PrunedTreeGenerator {
    
    Long seed=new Long(1);
    
    //list of new files
    public ArrayList<File> listPrunedAlignments=new ArrayList<>();
    public ArrayList<File> listPrunedTrees=new ArrayList<>();
    public File fileDtx=null;
    public File fileD2tx=null;
    
    //pruning percent
    double percentPruning=0.1; //10%
    
    
    //read generation: nomral distrib around mean R with sd (R/4)
    //and min length m
    int[] R={3,6};
    int m=2; //let's consider that we have at least 75bp reads

    //set which k/alpha are tested
    int k=5;
    int maxK=12;
    int kIncrement=1;
    double factor=1.0;
    double maxFactor=2.0;
    double factorIncrement=0.1;
    
    public static void main(String[] args) {
        
        try {
            
            ///////////////////////////////////////////////////////////////////
            //TEST ZONE
            
            //directory where all testes are done
            File workDir=new File("/media/ben/STOCK/DATA/viromeplacer/accu_tests");
            File treeFile=new File(workDir+File.separator+"RAxML_result.bv_refs_aln");
            File alignFile=new File(workDir+File.separator+"bv_refs_aln_stripped_99.5.fasta");
            
            
            
            //very small dataset for debugging
//            String HOME = System.getenv("HOME");
//            File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");
//            String treeString="(A:0.1,B:0.2,((C:0.1,D:0.2)Y:0.1,(E:0.1,F:0.2)Z:0.2)X:0.2)W:0.0;";
//            Fasta fA=new Fasta("A", "AAAAAT");
//            Fasta fB=new Fasta("B", "AAATAT");
//            Fasta fC=new Fasta("C", "GGAATA");
//            Fasta fD=new Fasta("D", "GGGATT");
//            Fasta fE=new Fasta("E", "CCATTC");
//            Fasta fF=new Fasta("F", "CCATGT"); 
//            ArrayList<Fasta> fastas=new ArrayList<>();
//            fastas.add(fA); fastas.add(fB); fastas.add(fC); fastas.add(fD);
//            fastas.add(fE); fastas.add(fF);
//            Alignment align=new Alignment(fastas);

            ////////////////////////////////////////
            //pplacer rRNA dataset
            
            //set directory associated to this tree
            workDir=new File(workDir.getAbsolutePath()+File.separator+"pplacer_16s");
            if (!workDir.exists()) {
                workDir.mkdir();
            }
            FASTAPointer fp=new FASTAPointer(alignFile, false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
            }
            Alignment align=new Alignment(fastas);
            System.out.println(align.describeAlignment(false));
            fp.closePointer();
            
            //read tree
            BufferedReader br=new BufferedReader(new FileReader(treeFile));
            String line=null;
            String treeString=null;
            while((line=br.readLine())!=null) {treeString=line;}
            PhyloTree tree = NewickReader.parseNewickTree2(treeString, true);
            tree.initIndexes();
            System.out.println("Original tree, # nodes: "+tree.getNodeCount());
            //tree.displayTree();

            
            //PREPARE ALL PRUNING EXPERIMENT FILES
            ///////////////////////////////////////////////////
            
            //launch
            PrunedTreeGenerator ptg=new PrunedTreeGenerator();

            ptg.generatePrunedTrees(workDir,tree,align);
            System.out.println("PRUNED TREE GENERATION DONE !");
            
            
            
            //LAUNCH DB_BUILD FOR ALL k/alpha COMBINATIONS
            ///////////////////////////////////////////////////
            ptg.buildDBs(workDir);
            System.out.println("BUILD OF ALL DB (different k/alpha) DONE !");
            
            
            
            
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
        tree.getNodeIdsByDFS().toArray(nodeIds);
        //shuffle their order
        shuffleArray(nodeIds);
        //define first x% as pruning experiments
        Integer[]prunedNodeIds=Arrays.copyOfRange(nodeIds, 0, new Double(percentPruning*nodeIds.length).intValue());
        System.out.println("pruning experiments: "+prunedNodeIds.length);
        System.out.println("prunedNodeIds: "+Arrays.toString(prunedNodeIds));
        //sort to make more comprehensive output matrices Dtx and D'tx
        Arrays.sort(prunedNodeIds);

        //log registering all skipped Nx, and the reason
        File skippedLog=new File(workDir+File.separator+"SKIPPED_Nx");
        BufferedWriter skippedLogBw=new BufferedWriter(new FileWriter(skippedLog,false));
        
        //prepare Dtx and D'tx matrices, build their headers
        fileDtx=new File(workDir+File.separator+"Dtx.csv");
        fileD2tx=new File(workDir+File.separator+"D2tx.csv");
        BufferedWriter brDtx=new BufferedWriter(new FileWriter(fileDtx));
        BufferedWriter brD2tx=new BufferedWriter(new FileWriter(fileD2tx));
        StringBuilder nodeDistMatrix=new StringBuilder();
        StringBuilder branchDistMatrix=new StringBuilder();
        int matrixSize=tree.getNodeCount();
        System.out.println("Dtx size (x,y): "+matrixSize+","+prunedNodeIds.length);
        nodeDistMatrix.append("nodeLabels;");
        for (int nodeId=0;nodeId<tree.getNodeCount();nodeId++)
            nodeDistMatrix.append(";"+tree.getById(nodeId).getLabel());
        nodeDistMatrix.append("\n");
        nodeDistMatrix.append(";nodeIds");
        for (int nodeId=0;nodeId<tree.getNodeCount();nodeId++)
            nodeDistMatrix.append(";"+nodeId);
        nodeDistMatrix.append("\n");
        branchDistMatrix.append(new String(nodeDistMatrix)); //simple contructor copy
        
        
            
        //launch pruning for each selected Nx
        for (int i = 0; i < prunedNodeIds.length; i++) {
            Integer nx_id = prunedNodeIds[i];
            System.out.println("--------------------------------------");
            System.out.println("copying tree and alignment before pruning...");
            PhyloNode rootCopy=tree.getRoot().copy();
            PhyloTree treeCopy=new PhyloTree(new PhyloTreeModel(rootCopy),tree.isRooted());;
            System.out.println("indexing tree ...");
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
            System.out.println("Starting pruning...");
            
            //current root Nx defining the pruned clade
            PhyloNode Nx= treeCopy.getById(nx_id);
            System.out.println("--------------------------------------");
            System.out.println("selected Nx: x="+i+"  node:"+Nx);
            if (Nx.isRoot()) {
                skippedLogBw.append("Nx="+i+"\t"+Nx.toString()+"\tIs root, so skipped.\n");
                System.out.println("SKIPPED: this node is root (id="+nx_id+", no pruning.");
                continue;
            }
                        
            
            
            
            ///////////////////////////////////////////////////////////
            //first, modify alignment and deleted subtree related to this Nx
            
            
            //if leaf, removes from multiple align
            //and export this leave to a it file
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
            //if this removal ends to a multiple alignment with only <3 leaves
            //we don't go further and pass to the next Nx.
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
            
            
            //prepare the corresponding output files
            File AxDir=new File(workDir+File.separator+"Ax");
            File TxDir=new File(workDir+File.separator+"Tx");
            File RxDir=new File(workDir+File.separator+"Rx");
            AxDir.mkdir();
            TxDir.mkdir();
            RxDir.mkdir();
            
            File Ax=new File(AxDir+File.separator+"A"+i+"_nx="+nx_id+"("+Nx.getLabel()+").align");
            File Tx=new File(TxDir+File.separator+"T"+i+"_nx="+nx_id+"("+Nx.getLabel()+").tree");
            listPrunedTrees.add(Tx);
            listPrunedAlignments.add(Ax);
            
            //prepare writers
            NewickWriter nw=new NewickWriter(Tx);            
            
            
            //write alignment to file
            alignCopy.writeAlignmentAsFasta(Ax);
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

            //memorize which edge is the best placement, this is a node, 
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
                treeCopy=new PhyloTree(new PhyloTreeModel(Np),tree.isRooted());
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
                //detach and get rid of Nx and Np
                Np.removeFromParent();
                Np=null;
            }
            
            
            //save the pruned tree
            System.out.println("Indexing pruned tree");
            treeCopy.initIndexes(); //usefull to remove former nodes from the maps
            //before being written
            System.out.println("Writing pruned tree newick");
            nw.writeNewickTree(treeCopy, true, true, false);  //no internal node names if PAML !
            ////////////////////////////////////////////////////////////////////
            //third, build Dtx and D'tx line corresponding to this pruning
            //note: for very big trees, should build that as an object
            //and save it by serialization ?
            nodeDistMatrix.append(Nx.getLabel()+";"+Nx.getId());
            branchDistMatrix.append(Nx.getLabel()+";"+Nx.getId());
            for (int n=0;n<prunedNodeIds.length;n++) {
                PhyloNode currentNode=treeCopy.getById(prunedNodeIds[n]);
                nodeDistMatrix.append(";");
                branchDistMatrix.append(";");
                if (currentNode==null) { //this node was pruned
                    nodeDistMatrix.append("-1");
                    branchDistMatrix.append("-1.0");
                    continue;
                }  else if (currentNode==Np_p) {
                    nodeDistMatrix.append("0");
                    branchDistMatrix.append("0.0");
                    //System.out.println("same node!");
                } else {
                    //retrieve path from thiw new edge to all other nodes
                    Path shortestPath = treeCopy.shortestPath(treeCopy.getRoot(), b_new,currentNode);
                    //System.out.println("path: "+shortestPath);
                    nodeDistMatrix.append(shortestPath.nodeDistance);
                    //for the branch ditance, we need to add the distance that
                    //was separating Np and Np'' if the path startt to go up
                    //between Np and Np' if the path start to go down
                    //
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

            for (int j = 0; j < R.length; j++) {
                //prepare the corresponding output files
                File Rxj=new File(RxDir+File.separator+"R"+i+"_nx="+nx_id+"("+Nx.getLabel()+")_r"+R[j]+".fasta");
                FileWriter fw=new FileWriter(Rxj);
                int r = R[j];
                for (Iterator<Fasta> it = leavesRemoved.iterator(); it.hasNext();) {
                    Fasta next = it.next();
                    int fasta_length=next.getSequence().length();
                    //select a normally distibution read length
                    ThreadLocalRandom randomGenerator = ThreadLocalRandom.current();
                    int v_length = new Double(r+Math.floor(randomGenerator.nextGaussian()*r/4)).intValue();
                    //System.out.println(v_length);
                    if (v_length>=m && v_length<fasta_length) {
                        //select an uniformely distibution position
                        int p=randomGenerator.nextInt(0, fasta_length-v_length);
                        String read=next.getSequence().substring(p, p+v_length);
                        Fasta n=new Fasta(next.getHeader()+"_r"+R[j]+"_"+p+":"+(p+v_length-1), read);
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
        
        
        //closes everything
        brDtx.append(nodeDistMatrix);
        brD2tx.append(branchDistMatrix);
        brDtx.close();
        brD2tx.close();
        skippedLogBw.close();
                
    }
    
    
    
    private void buildDBs(File workDir) {

        ArrayList<File> workDirs=new ArrayList<>();
        for (int i = 0; i < listPrunedAlignments.size(); i++) {
            File a=listPrunedAlignments.get(i);
            String experiment=a.getName().split("\\.")[0];
            System.out.println("LAUCHING ALL DB_BUILDs FOR EXPERIMENT: "+experiment);
            File DBx=new File(workDir.getAbsolutePath()+File.separator+"Dx");
            DBx.mkdir();
            File expPath=new File(DBx+File.separator+experiment);
            System.out.println("Experiment work dir: "+expPath.getAbsolutePath());
            workDirs.add(expPath);
            if (workDir.canWrite())
                expPath.mkdir();
            else {
                System.out.println("Cannot write in dir: "+expPath.getAbsolutePath());
                System.exit(1);
            }
            NumberFormat nf = NumberFormat.getNumberInstance();
            nf.setMaximumFractionDigits(1);

            for (int current_k = k; current_k < maxK+kIncrement; current_k+=kIncrement) {
                for (double current_alpha = factor; current_alpha < maxFactor+factorIncrement; current_alpha+=factorIncrement) {

                    Infos.println("#######################################################");
                    Infos.println("# NEW DB_BUILD LAUNCH; Config: k="+current_k+" factor="+current_alpha);
                    Infos.println("#######################################################");

                    //make subdirectory corresponding to this k/alpha combination
                    File combDir=new File(expPath.getAbsolutePath()+File.separator+"k"+current_k+"_a"+nf.format(current_alpha));
                    boolean dirOK=combDir.mkdir();
                    System.out.println("Current work_dir: "+combDir.getAbsolutePath());

                    ////////////////////////////////////////////////////////
                    //need to change the approach !!!
                    //-add an option to viromeplacer to use an already existing
                    // ancestral reconstruction from a directory
                    //-build the baseml.ctl in this diectory, using the 
                    //current script
                    //-then generate all the qsub commands necessary to launch
                    //many qsub jobs, one per paml AR. With vaginal dataset
                    //could do all in just 1 hour.
                    //
                    /////////////////////////////////////////////////////////
                    
                    
                    if (dirOK) {
                        List<String> com=new ArrayList<>();
                        //com.add("/usr/bin/java");
                        //com.add("-jar");
                        //com.add("/home/benclaff/NetBeansProjects/viromeplacer/dist/ViromePlacer.jar");
                        //com.add("/media/ben/STOCK/SOURCES/NetBeansProjects/ViromePlacer/dist/ViromePlacer.jar");
                        com.add("-m");
                        com.add("b");
                        com.add("-a");
                        com.add(new Float(current_alpha).toString());
                        com.add("-k");
                        com.add(new Integer(current_k).toString());
                        com.add("-t");
                        com.add(listPrunedTrees.get(i).getAbsolutePath());
                        com.add("-i");
                        com.add(listPrunedAlignments.get(i).getAbsolutePath());
                        com.add("-w");
                        com.add(combDir.getAbsolutePath());
                        com.add("-v");
                        com.add("0");
                        //executeProcess(com, expPath);
                        String[] arguments = new String[com.size()];
                        arguments = com.toArray(arguments);
                        Main_v2.main(arguments);

                        //DBGeneration(fw,current_k,(float)current_alpha);
                        System.gc();
                    } else {
                        System.out.println("AR reconstruction directory could not be created: "+combDir.getAbsolutePath());
                        System.exit(1);
                    }
                    
                }

            }



        }
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
    
}
