
import alignement.Alignment;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.attribute.PosixFilePermission;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
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
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");

    //list of new files
    public List<File> prunedAlignmentsFiles=new ArrayList<>(); //list of pruned alignments
    public List<File> prunedTreesFiles=new ArrayList<>(); //list of pruned Trees
   
    //where is EPA 
    File RAxMLBinary=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/RAxML-8.2.9/raxmlHPC-PTHREADS-SSE3");

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
        
        FileWriter fw=null;
        try {
            System.out.println("ARGS: workDir RAxMLBinary");
            
            //launch
            RAXMLExperiment exp=new RAXMLExperiment();
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                exp.workDir=new File(args[0]);
                exp.RAxMLBinary=new File(args[1]);
                System.out.println("workDir: "+exp.workDir);
                System.out.println("RAxMLBinary: "+exp.RAxMLBinary);
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
            }   //build directory /RAXMLx
            File EPAxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"EPAx");
            EPAxDir.mkdir();
            File HMMxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"HMMx");
            //list of raxml commands
            StringBuilder sbQsubCommands=new StringBuilder();
            //build all Ax directories in /HMMx and for each set of hmm-aligned
            //reads (/HMMx/Ax_nx_la/Rx_nx_la_r), build the list of commands
            for (int i=0;i<exp.prunedAlignmentsFiles.size();i++) {
                String experimentLabel=exp.prunedAlignmentsFiles.get(i).getName().split("\\.")[0];
                System.out.println("EPA commands for "+experimentLabel);
                File EPAxAxDir=new File(EPAxDir.getAbsolutePath()+File.separator+experimentLabel);
                EPAxAxDir.mkdir();
                File HMMxAxDir=new File(HMMxDir.getAbsolutePath()+File.separator+experimentLabel);
                System.out.println("Searching alignments in "+HMMxAxDir.getAbsolutePath());
                File[] listFiles = HMMxAxDir.listFiles((File dir, String name) -> {
                    if (name.endsWith(".aln.fasta")) {
                        return true;
                    } else {
                        return false;
                    }
                });
                for (int j = 0; j < listFiles.length; j++) {
                    File f = listFiles[j];
                    String readAlignLabel=f.getName().split("\\.")[0];
                    //build the placement commands
                    //example:  RAxMLBinary -f v -G 0.1 -m GTRCAT -n EPA -s ../concat_basic_queries.fasta -t ../RAxML_bestTree.basic_tree
                    StringBuilder sbRAxMLCommand=new StringBuilder();
                    sbRAxMLCommand.append(  exp.RAxMLBinary.getAbsolutePath()+" " +
                            "-f v -G 0.1 -m GTRCAT " +
                            "-n "+readAlignLabel+" " +
                            "-s "+f.getAbsolutePath()+" " +
                            "-t "+exp.prunedTreesFiles.get(i).getAbsolutePath() +
                            "\n"
                    );
                    //do its qsub counterpart
                    sbQsubCommands.append("echo \""+sbRAxMLCommand.toString()+"\" |  qsub -N EPAx_"+readAlignLabel+
                            " -wd "+EPAxAxDir.getAbsolutePath()+"\n");
                }
            }
            File qsubEPACommands=new File(EPAxDir.getAbsolutePath()+File.separator+"qsub_epa_commands");
            fw = new FileWriter(qsubEPACommands);
            fw.append(sbQsubCommands.toString());
            fw.close();
            Files.setPosixFilePermissions(qsubEPACommands.toPath(), exp.perms);
            
        } catch (IOException ex) {
            Logger.getLogger(RAXMLExperiment.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(RAXMLExperiment.class.getName()).log(Level.SEVERE, null, ex);
            }  
        }
            

    }  
}
