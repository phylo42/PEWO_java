
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.NewickReader;
import tree.NewickWriter;
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
public class GenerateAllTrifurcationsOfUnrootedTree {
    
    public static void main(String[] args) {
        
        try {
            
            //read original file
            File treeFile = new File("/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/pplacer_16S_rRNA/RAxML_result.bv_refs_aln");
            String s=null;
            BufferedReader brInputNewick = Files.newBufferedReader(treeFile.toPath());
            //file should contain a single newick
            while ( (s=brInputNewick.readLine())!=null ) {break;}
            PhyloTree tree = NewickReader.parseNewickTree2(s, false, false);
            tree.initIndexes();
            
            //only allow unrooted trees
            if (tree.isRooted()) {
                System.out.println("Tree is rooted, abort.");
                System.exit(1);
            }
            
            ArrayList<Integer> nodeIdsByDFS = tree.getNodeIdsByDFS();
            int numberTrifurcationsExpected=tree.getLeavesCount()-3;
            
            System.out.println("Will generate "+numberTrifurcationsExpected+" trifurcations. ");

            int counter=0;
            for (int i = 0; i < nodeIdsByDFS.size(); i++) {
                Integer nodeId = nodeIdsByDFS.get(i);
                
                if (!tree.getById(nodeId).isLeaf()) {
                    
                    //do not reroot on original root
                    if (tree.getById(nodeId).isRoot()) {
                        System.out.println("Skipping node of original root.");
                        continue;
                    }
                    
                    if ((counter%100)==0) {
                        System.out.println(i+"/"+numberTrifurcationsExpected);
                    }
                    
                    //copy tree
                    PhyloNode rootCopy=tree.getRoot().copy();
                    PhyloTree treeCopy=new PhyloTree(new PhyloTreeModel(rootCopy),tree.isRooted(), false);
                    treeCopy.initIndexes();
                    //treeCopy.displayTree();
                    //if (true) break;
                    
                    //reroot on nodeId
                    PhyloNode node=treeCopy.getById(nodeId);
                    System.out.println("Rerooting on node: "+node.toString());
                    treeCopy.rerootTree(treeCopy.getById(nodeId), false);
                                        
                    
                    //output as newick
                    File outname=new File("/home/ben/Dropbox/viromeplacer/test_datasets/trifurcations_test/reroot_"+counter+".tree");
                    BufferedWriter bwOut = Files.newBufferedWriter(outname.toPath());
                    new NewickWriter(bwOut).writeNewickTree(treeCopy, true, true, false, false);
                    bwOut.close();

                    //discard copy
                    treeCopy=null; 
                    
                    counter++;
                }
                
            }
            
            System.out.println("# trifurcations generated: "+counter);
            
            brInputNewick.close();
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(GenerateAllTrifurcationsOfUnrootedTree.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
}
