
import alignement.Alignment;
import inputs.Fasta;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
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
            
            
            PrunedTreeGenerator ptg=new PrunedTreeGenerator();
            ptg.generatePrunedTrees(workDir,tree,align);
            
            Thread.sleep(10000);
            
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
            //some checkup about the original  copy
            System.out.println("Tree; #nodes="+treeCopy.getNodeCount());
            System.out.println("Tree; rooted="+treeCopy.isRooted());
            System.out.println("Tree; #leaves="+treeCopy.getLeavesCount());
            Enumeration depthFirstEnumeration = treeCopy.getRoot().depthFirstEnumeration();
            while(depthFirstEnumeration.hasMoreElements()) {
                System.out.println("Tree;  nodes="+depthFirstEnumeration.nextElement());
            }
            System.out.println("Starting pruning...");
            //current root Nx defininf the pruned clade
            PhyloNode Nx= treeCopy.getById(nx_id);
            System.out.println("--------------------------------------");
            System.out.println("selected Nx: "+Nx+" (id="+nx_id+")");
            
            
            //TO DO : constructor copy of the alignment, like the tree
            
            
            
            ///////////////////////////////////////////////////////////
            //first, modify alignment and deleted subtree related to this Nx
            
            //prepare the corresponding output files
            File Ax=new File(workDir+File.separator+"A"+i+"_nx="+nx_id+"("+Nx.getLabel()+").align");
            File Tx=new File(workDir+File.separator+"T"+i+"_nx="+nx_id+"("+Nx.getLabel()+").tree");
            File Dtx=new File(workDir+File.separator+"Dt"+i+"_nx="+nx_id+"("+Nx.getLabel()+").csv");
            File D2tx=new File(workDir+File.separator+"D2t"+i+"_nx="+nx_id+"("+Nx.getLabel()+").csv");
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
                align.removeSequence(Nx.getLabel());
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
                        align.removeSequence(nextNx.getLabel());
                    } 
                }
                //delete Nx children (this does not delete the sub nodes from memory !!!)
                nextNx.removeAllChildren();
            }
            //write alignment to file
            align.writeAlignmentAsFasta(Ax);
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
            //save the pruned tree
            System.out.println("Indexing pruned tree");
            treeCopy.initIndexes(); //usefull to remove former nodes from the maps
            //before being written
            System.out.println("Writing pruning tree newick");
            nw.writeNewickTree(treeCopy, true, true, false);
            

            
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
    
    
    
}
