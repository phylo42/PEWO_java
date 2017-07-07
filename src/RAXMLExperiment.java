
import alignement.Alignment;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import tree.NewickReader;
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
public class RAXMLExperiment {
    Long seed=new Long(1);
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");

    //list of new files
    public List<File> prunedAlignmentsFiles=new ArrayList<>(); //list of pruned alignments
    public List<File> prunedTreesFiles=new ArrayList<>(); //list of pruned Trees
    public File fileDtx=null;
    public File fileD2tx=null;
    
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

    //where is EPA 
    File RAxMLBinary=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/RAxML-8.2.9/raxmlHPC-PTHREADS-SSE3");
    
    
    public static void main(String[] args) {
        
        System.out.println("ARGS: workDir");
        
            //launch
            RAXMLExperiment exp=new RAXMLExperiment();
            
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                exp.workDir=new File(args[0]);
                System.out.println("workDir: "+exp.workDir);
            }
            
            //load alignments trees from Ax/Tx directories
            File Ax=new File(exp.workDir.getAbsolutePath()+File.separator+"Ax");
            File Tx=new File(exp.workDir.getAbsolutePath()+File.separator+"Tx");
            exp.prunedAlignmentsFiles=Arrays.stream(Ax.listFiles()).sorted().collect(Collectors.toList());
            exp.prunedTreesFiles=Arrays.stream(Tx.listFiles()).sorted().collect(Collectors.toList());
            
            System.out.println(exp.prunedAlignmentsFiles);
            System.out.println(exp.prunedTreesFiles);

            if (    exp.prunedAlignmentsFiles.size()<1 ||
                    exp.prunedTreesFiles.size()<1 ||
                    exp.prunedAlignmentsFiles.size()!=exp.prunedTreesFiles.size()) {
                System.out.println("Unexpected Ax/Tx files");
                System.exit(1); 
            }
            
            //build directory /RAXMLx
            File EPAxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"EPAx");
            EPAxDir.mkdir();
            //build all pruning directories in /RAXMLx and build the list of commands
            for (int i=0;i<exp.prunedAlignmentsFiles.size();i++) {
                String experimentLabel=exp.prunedAlignmentsFiles.get(i).getName().split("\\.")[0];
                File expDir=new File(EPAxDir.getAbsolutePath()+File.separator+experimentLabel);
                expDir.mkdir();
                
                //build the hmmalign command
                
                //#build hmm profile from multiple alignment (.aln file in fasta)
                //hmmbuild basic.hmm basic/basic.aln 
                //../../../hmmer-3.1b2/binaries/hmmbuild A0_nx4_la.hmm ../../Ax/A0_nx4_la.align
                
                //#align reads using Hmm profile
                //hmmalign -o basic.sto --mapali basic/basic.aln basic.hmm ../queries.fasta
                //(working? --outformat FASTA ; no! we should use PSIBLAST output format and rapidly convert it to fasta --outformat PSIBLAST)
                //../../../hmmer-3.1b2/binaries/hmmalign --outformat PSIBLAST -o ./A0_nx4_la.reads.align --mapali ../../Ax/A0_nx4_la.align A0_nx4_la.hmm ../../Rx/R0_nx4_la_r300.fasta

                
                //build the placement command
                //  RAxMLBinary -f v -G 0.1 -m GTRCAT -n EPA -s ../concat_basic_queries.fasta -t ../RAxML_bestTree.basic_tree 
                
            }
            
            
            
            
            
    }  
}
