
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.attribute.PosixFilePermission;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class PPLACERExperiment {
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");

    //list of new files
    public List<File> prunedAlignmentsFiles=new ArrayList<>(); //list of pruned alignments
    public List<File> prunedTreesFiles=new ArrayList<>(); //list of pruned Trees
   
    //where where are taxit and pplacer  
    File TAXITBinary=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/pplacer-Linux-v1.1.alpha18-2-gcb55169/taxtastic/taxit");
    File PPLACERBinary=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/pplacer-Linux-v1.1.alpha18-2-gcb55169/pplacer");
    
    //stats files, required by taxit
    File statFile=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/RAxML_info.6_leaves_test_set");


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
            System.out.println("ARGS: workDir PPLACERBinary TAXITBinary statisticFile(ex:RAxML_info_xxx)");
            
            //launch
            PPLACERExperiment exp=new PPLACERExperiment();
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                exp.workDir=new File(args[0]);
                exp.PPLACERBinary=new File(args[1]);
                exp.TAXITBinary=new File(args[2]);
                exp.statFile=new File(args[3]);
                System.out.println("workDir: "+exp.workDir);
                System.out.println("PPLACERBinary: "+exp.PPLACERBinary);
                System.out.println("TAXITBinary: "+exp.TAXITBinary);
                System.out.println("statFile: "+exp.statFile);
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
            File PPLxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"PPLx");
            PPLxDir.mkdir();
            File HMMxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"HMMx");
            //list of raxml commands
            StringBuilder sbQsubCommands=new StringBuilder();
            //build all Ax directories in /HMMx and for each set of hmm-aligned
            //reads (/HMMx/Ax_nx_la/Rx_nx_la_r), build the list of commands
            for (int i=0;i<exp.prunedAlignmentsFiles.size();i++) {
                String experimentLabel=exp.prunedAlignmentsFiles.get(i).getName().split("\\.")[0];
                System.out.println("PPL commands for "+experimentLabel);
                File PPLxAxDir=new File(PPLxDir.getAbsolutePath()+File.separator+experimentLabel);
                PPLxAxDir.mkdir();
                File HMMxAxDir=new File(HMMxDir.getAbsolutePath()+File.separator+experimentLabel);
                System.out.println("Searching alignments in "+HMMxAxDir.getAbsolutePath());
                File[] listFiles = HMMxAxDir.listFiles((File dir, String name) -> {
                    if (name.endsWith(".aln.fasta")) {
                        return true;
                    } else {
                        return false;
                    }
                });
                //build the placement commands
                //will output both the taxit and pplacer commands in a signle script file
                File script=new File(PPLxAxDir.getAbsolutePath()+File.separator+"pplacer_pipeline.sh");
                FileWriter fwScript=new FileWriter(script);
                StringBuilder sbPPlacerCommands=new StringBuilder();
                //create pplacer database from this pruned alignment/tree
                //example:  taxit create -P ./basic -l basic -f ../basic.aln -t ../RAxML_bestTree.basic_tree -s ../RAxML_info.basic_tree 
                sbPPlacerCommands.append(  exp.TAXITBinary.getAbsolutePath()+" create" +
                                        " -P " +PPLxAxDir.getAbsolutePath()+File.separator+"refpkg" +
                                        " -l locus" +
                                        " -f " + exp.prunedAlignmentsFiles.get(i).getAbsolutePath() +
                                        " -t " + exp.prunedTreesFiles.get(i).getAbsolutePath() +
                                        " -s " + exp.statFile.getAbsolutePath() +
                                        "\n"
                );
                //for each read file, do a placement using the pplacer database
                for (int j = 0; j < listFiles.length; j++) {
                    File f = listFiles[j];
                    String readAlignLabel=f.getName().split("\\.")[0];

                    //example:  pplacer --verbosity 2 -c basic basic.sto 
                    sbPPlacerCommands.append(  exp.PPLACERBinary.getAbsolutePath() +
                                            " --verbosity 2" +
                                            " -c " +PPLxAxDir.getAbsolutePath()+File.separator+"refpkg" +
                                            " "+f.getAbsolutePath() +
                                            "\n"
                    );
                }
                //put these commands in a file
                fwScript.append(sbPPlacerCommands);
                fwScript.close();
                Files.setPosixFilePermissions(script.toPath(), exp.perms);
                //execute script through qsub
                sbQsubCommands.append(  "echo \""+script.getAbsolutePath()+"\" |" +
                                        " qsub -N PPLx_"+experimentLabel+
                                        " -wd "+PPLxAxDir.getAbsolutePath()+"\n");
            }
            File qsubEPACommands=new File(PPLxDir.getAbsolutePath()+File.separator+"qsub_pplacer_commands");
            fw = new FileWriter(qsubEPACommands);
            fw.append(sbQsubCommands.toString());
            Files.setPosixFilePermissions(qsubEPACommands.toPath(), exp.perms);
            
        } catch (IOException ex) {
            Logger.getLogger(PPLACERExperiment.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(PPLACERExperiment.class.getName()).log(Level.SEVERE, null, ex);
            }  
        }
            

    }   
}
