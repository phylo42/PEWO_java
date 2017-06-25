
import alignement.Alignment;
import etc.Infos;
import inputs.Fasta;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import main_v2.ARProcessLauncher;
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
    File workDir=null;
    
    //list of new files
    public ArrayList<File> listPrunedAlignments=new ArrayList<>();
    public ArrayList<File> listPrunedTrees=new ArrayList<>();
    public ArrayList<File> listDtx=new ArrayList<>();
    public ArrayList<File> listD2Tx=new ArrayList<>();  
    
    public static void main(String[] args) {
        
        try {
            
            ///////////////////////////////////////////////////////////////////
            //TEST ZONE, forces arguments
            String HOME = System.getenv("HOME");
            File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests");
            
            //load a test tree
            String treeString="(A:0.1,B:0.2,((C:0.1,D:0.2)Y:0.1,(E:0.1,F:0.2)Z:0.2)X:0.2)W:0.0;";
            
            Fasta fA=new Fasta("A", "AAAAA");
            Fasta fB=new Fasta("B", "AAATA");
            Fasta fC=new Fasta("C", "GGAAT");
            Fasta fD=new Fasta("D", "GGGAT");
            Fasta fE=new Fasta("E", "CCATT");
            Fasta fF=new Fasta("F", "CCATG"); 
            ArrayList<Fasta> fastas=new ArrayList<>();
            fastas.add(fA); fastas.add(fB); fastas.add(fC); fastas.add(fD);
            fastas.add(fE); fastas.add(fF);
            Alignment align=new Alignment(fastas);
            
            PhyloTree tree = NewickReader.parseNewickTree2(treeString, true);
            tree.initIndexes();
            
            //tree.displayTree();
            ////////////////////////////////////////////////////////
            
          
            
            PrunedTreeGenerator ptg=new PrunedTreeGenerator();
            ptg.generatePrunedTrees(workDir,tree,align);
            System.out.println("PRUNED TREE GENERATION DONE !");
            
            //launching DB_BUILD on all these files
            ArrayList<File> workDirs=new ArrayList<>();
            for (int i = 0; i < ptg.listPrunedAlignments.size(); i++) {
                File a=ptg.listPrunedAlignments.get(i);
                String experiment=a.getName().split("\\.")[0];
                System.out.println("LAUCHING DB_BUILD FOR EXPERIMENT: "+experiment);
                File expPath=new File(workDir.getAbsolutePath()+File.separator+experiment);
                System.out.println("Experiment work dir: "+expPath.getAbsolutePath());
                workDirs.add(expPath);
                if (workDir.canWrite())
                    expPath.mkdir();
                else {
                    System.out.println("Cannot write in dir: "+expPath.getAbsolutePath());
                    System.exit(1);
                }
                
                float alpha=1.0f;
                int k=3;
                
                List<String> com=new ArrayList<>();
                com.add("/usr/bin/java");
                com.add("-jar");
                com.add("/home/benclaff/NetBeansProjects/viromeplacer/dist/ViromePlacer.jar");
                com.add("-m");
                com.add("b");
                com.add("-a");
                com.add(new Float(alpha).toString());
                com.add("-k");
                com.add(new Integer(k).toString());
                com.add("-t");
                com.add(ptg.listPrunedTrees.get(i).getAbsolutePath());
                com.add("-i");
                com.add(ptg.listPrunedAlignments.get(i).getAbsolutePath());
                com.add("-w");
                com.add(expPath.getAbsolutePath());
                com.add("-v");
                com.add("1");
                               
                
                executeProcess(com, expPath);
                
                
                
            }
            
            
            
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
        //define 10 first %
        double percent=1;
        Integer[]prunedNodeIds=Arrays.copyOfRange(nodeIds, 0, new Double(percent*nodeIds.length).intValue());
        System.out.println("prunedNodeIds: "+Arrays.toString(prunedNodeIds));
        
        //launch pruning for each selected Nx
        for (int i = 0; i < prunedNodeIds.length; i++) {
            Integer nx_id = prunedNodeIds[i];
            System.out.println("--------------------------------------");
            System.out.println("copying tree before pruning...");
            PhyloNode rootCopy=tree.getRoot().copy();
            PhyloTree treeCopy=new PhyloTree(new PhyloTreeModel(rootCopy),tree.isRooted());;
            System.out.println("indexing tree ...");
            treeCopy.initIndexes();
            System.out.println("copy done");
//            //some checkup about the original  copy
//            System.out.println("Tree; #nodes="+treeCopy.getNodeCount());
//            System.out.println("Tree; rooted="+treeCopy.isRooted());
//            System.out.println("Tree; #leaves="+treeCopy.getLeavesCount());
//            Enumeration depthFirstEnumeration = treeCopy.getRoot().depthFirstEnumeration();
//            while(depthFirstEnumeration.hasMoreElements()) {
//                System.out.println("Tree;  nodes="+depthFirstEnumeration.nextElement());
//            }
            //copying alignment
            System.out.println("Copying alignment");
            Alignment alignCopy=align.copy();
            
            System.out.println("Starting pruning...");
            //current root Nx defininf the pruned clade
            PhyloNode Nx= treeCopy.getById(nx_id);
            System.out.println("--------------------------------------");
            System.out.println("selected Nx: "+Nx+" (id="+nx_id+")");
            if (Nx.isRoot()) {
                System.out.println("SKIPPED: this node is root (id="+nx_id+", no pruning.");
                continue;
            }
                        
            
            
            
            ///////////////////////////////////////////////////////////
            //first, modify alignment and deleted subtree related to this Nx
            
            //prepare the corresponding output files
            File Ax=new File(workDir+File.separator+"A"+i+"_nx="+nx_id+"("+Nx.getLabel()+").align");
            File Tx=new File(workDir+File.separator+"T"+i+"_nx="+nx_id+"("+Nx.getLabel()+").tree");
            File Dtx=new File(workDir+File.separator+"Dt"+i+"_nx="+nx_id+"("+Nx.getLabel()+").csv");
            File D2tx=new File(workDir+File.separator+"D2t"+i+"_nx="+nx_id+"("+Nx.getLabel()+").csv");
            listPrunedTrees.add(Tx);
            listPrunedAlignments.add(Ax);
            listDtx.add(Dtx);
            listD2Tx.add(D2tx);
            //prepare writers
            BufferedWriter brDtx=new BufferedWriter(new FileWriter(Dtx));
            BufferedWriter brD2tx=new BufferedWriter(new FileWriter(D2tx));
            NewickWriter nw=new NewickWriter(Tx);            
            
            //if leaf, removes from multiple align
            //and export this leave to a fasta file
            ArrayList<String> leavesRemoved=new ArrayList<>();
            if (Nx.isLeaf()) {
                System.out.println("nx is leaf !");
                leavesRemoved.add(Nx.getLabel());
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
                    System.out.println("DFS son removed "+nextNx);
                    //if leaf, removes and export
                    if (nextNx.isLeaf()) {
                        leavesRemoved.add(nextNx.getLabel());
                        alignCopy.removeSequence(nextNx.getLabel());
                    } 
                }
                //delete Nx children (this does not delete the sub nodes from memory !!!)
                nextNx.removeAllChildren();
            }
            //write alignment to file
            alignCopy.writeAlignmentAsFasta(Ax);
            System.out.println("LeavesRemoved: "+leavesRemoved);
            
            ////////////////////////////////////////////////////////////////////
            //second, built Tx
            //remove selected Nx itself and link its sister to parent
            PhyloNode Np=(PhyloNode)Nx.getParent(); //Np
            System.out.println("Np:"+Np);
            float Np_bl=Np.getBranchLengthToAncestor(); //bl of Np to Np''
            PhyloNode Np_pp=(PhyloNode)Np.getParent(); //Np''
            System.out.println("Np_pp "+Np_pp);
            PhyloNode Np_p=(PhyloNode)Np.getChildBefore(Nx); //Np' if on left
            if (Np_p==null) //if not on left, will be on right
                Np_p=(PhyloNode)Np.getChildAfter(Nx); //Np' if on right
            System.out.println("Np_p "+Np_p);
            float Np_p_bl=Np_p.getBranchLengthToAncestor(); //bl of Np' to Np
            
            //when Nx is a root's son, Np_pp (parent of root) doesn't exists
            //basically, in this case we just remove the subtree of Nx and
            // the root; Np' becomes the root
            if (Np_pp==null) { 
                Np.removeFromParent();//does nothing in fact
                Np_p.removeFromParent();
                Nx.removeFromParent();
                Np_p.setBranchLengthToAncestor(0.0f);
                //rebuilt phylotree using this new root
                treeCopy=new PhyloTree(new PhyloTreeModel(Np_p),tree.isRooted());
                
            //in all other case, the tree is disconnected and reconnected
            } else{ 
                //disconnect all
                Np.removeFromParent();
                Np_p.removeFromParent();            
                Nx.removeFromParent();
                //connect Np' tp Np'' and update bl
                Np_pp.add(Np_p);
                Np_p.setBranchLengthToAncestor(Np_bl+Np_p_bl);
                //detach and get rid of Nx and Np
                Np.removeFromParent();
                Np=null;
                Nx.removeFromParent();
                Nx=null;
            }
            
            
            //save the pruned tree
            System.out.println("Indexing pruned tree");
            treeCopy.initIndexes(); //usefull to remove former nodes from the maps
            //before being written
            System.out.println("Writing pruning tree newick");
            nw.writeNewickTree(treeCopy, true, true, false);  //no internal node names if PAML !
            

            
            ////////////////////////////////////////////////////////////////////
            //third, build Dtx and D'tx in a single go
            //note: for very big trees, should build that as an object
            //and save it by serialization
            StringBuilder nodeDistMatrix=new StringBuilder();
            StringBuilder branchDistMatrix=new StringBuilder();
            int matrixSize=treeCopy.getNodeCount();
            ArrayList<Integer> nodeIdsByDFS = treeCopy.getNodeIdsByDFS();
            System.out.println("Dtx Nodes:"+treeCopy.getNodeIdsByDFS());
            System.out.println("Dtx size (x,y):"+matrixSize);
            
            //header of both matrices
            nodeDistMatrix.append("nodeLabels;");
            for (int nodeId:nodeIdsByDFS)
                nodeDistMatrix.append(";"+treeCopy.getById(nodeId).getLabel());
            nodeDistMatrix.append("\n");
            nodeDistMatrix.append(";nodeIds");
            for (int nodeId:nodeIdsByDFS)
                nodeDistMatrix.append(";"+nodeId);
            nodeDistMatrix.append("\n");
            branchDistMatrix.append(new String(nodeDistMatrix)); //simple contructor copy
            
            //distances
            for (int y = 0; y < matrixSize; y++) {
                int yId=nodeIdsByDFS.get(y);
                PhyloNode nodeY=treeCopy.getById(yId);
                nodeDistMatrix.append(nodeY.getLabel()+";"+yId);
                branchDistMatrix.append(nodeY.getLabel()+";"+yId);
                for (int x = 0; x < matrixSize; x++) {
                    PhyloNode nodeX=treeCopy.getById(nodeIdsByDFS.get(x));
                    nodeDistMatrix.append(";");
                    branchDistMatrix.append(";");
                    if (nodeX==nodeY) {
                        nodeDistMatrix.append("0");
                        branchDistMatrix.append("0.0");
                    }
                    else {
                        Path shortestPath = treeCopy.shortestPath(treeCopy.getRoot(), nodeX,nodeY);
                        nodeDistMatrix.append(shortestPath.nodeDistance);
                        branchDistMatrix.append(shortestPath.branchDistance);
                    }
                }
                nodeDistMatrix.append("\n");
                branchDistMatrix.append("\n");
            }
            //System.out.println("----------------");
            //System.out.println(nodeDistMatrix.toString());
            //System.out.println("----------------");
            //System.out.println(branchDistMatrix.toString());
            //System.out.println("----------------");
            brDtx.append(nodeDistMatrix);
            brD2tx.append(branchDistMatrix);
          

            
            //debug
            //tree.displayTree();
            
            //closes everything
            brDtx.close();
            brD2tx.close();
            nw.close();
            treeCopy=null;
            
            
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
        FileOutputStream STDOUTOutputStream=new FileOutputStream(new File(workDir.getAbsolutePath()+File.separator+"sdtout.txt"));
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
        t.setDaemon(true);
        t.start();
    }
    
}
