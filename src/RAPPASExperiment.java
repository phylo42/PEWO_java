
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
public class RAPPASExperiment {
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");

    //list of new files
    public List<File> prunedAlignmentsFiles=new ArrayList<>(); //list of pruned alignments
    public List<File> prunedTreesFiles=new ArrayList<>(); //list of pruned Trees
   
    //where is EPA 
    File RAPPAJar=new File(HOME+"/home/ben/Dropbox/viromeplacer/test_datasets/software/rappas/ViromePlacer.jar");

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
            System.out.println("ARGS: workDir RAPPASJar");
            
            //launch
            RAPPASExperiment exp=new RAPPASExperiment();
            
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                exp.workDir=new File(args[0]);
                exp.RAPPAJar=new File(args[1]);
                System.out.println("workDir: "+exp.workDir);
                System.out.println("RAPPASJar: "+exp.RAPPAJar);
            }  
            

            //load alignments trees from AxDir/Tx directories
            File AxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"Ax");
            File TxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"Tx");
            exp.prunedAlignmentsFiles=Arrays.stream(AxDir.listFiles()).sorted().collect(Collectors.toList());
            exp.prunedTreesFiles=Arrays.stream(TxDir.listFiles()).sorted().collect(Collectors.toList());
            System.out.println(exp.prunedAlignmentsFiles);
            System.out.println(exp.prunedTreesFiles);
            if (    exp.prunedAlignmentsFiles.size()<1 ||
                    exp.prunedTreesFiles.size()<1 ||
                    exp.prunedAlignmentsFiles.size()!=exp.prunedTreesFiles.size()) {
                System.out.println("Unexpected Ax/Tx files");
            }
            
            //load list of k/alpha combinations
            String firstExperimentLabel=exp.prunedAlignmentsFiles.get(0).getName().split("\\.")[0];
            File DxDir=new File(exp.workDir.getCanonicalPath()+File.separator+"Dx");
            File DxA0Dir=new File(exp.workDir.getCanonicalPath()+File.separator+"Dx"+File.separator+firstExperimentLabel);
            File[] listFiles = DxA0Dir.listFiles();
            ArrayList<File> kAlphaDirs=new ArrayList<>();
            ArrayList<String> kAlphaDirsNames=new ArrayList<>();
            for (int i = 0; i < listFiles.length; i++) {
                File kAlphaDir = listFiles[i];
                if (kAlphaDir.getName().startsWith("k")) {
                    kAlphaDirs.add(kAlphaDir);
                    kAlphaDirsNames.add(kAlphaDir.getName());
                }
            }
            System.out.println(kAlphaDirsNames);
            
            ////////////////////////////////////////////////////////////////
            //first build build_db commands
            
            //list of rappas db_build commands
            StringBuilder sbQsubDBBUILDCommands=new StringBuilder();
            File qsub_dbbuild_commands=new File(DxDir.getAbsolutePath()+File.separator+"qsub_dbbuild_commands");
            //for all AxDir, each combination of k/alpha
            for (int i = 0; i < exp.prunedAlignmentsFiles.size(); i++) {
                File Ax = exp.prunedAlignmentsFiles.get(i);
                File Tx = exp.prunedTreesFiles.get(i);
                String experimentLabel=Ax.getName().split("\\.")[0];
                for (int j = 0; j < kAlphaDirsNames.size(); j++) {
                    String  DxAxKAlphaDirName= kAlphaDirsNames.get(j);
                    File  DxAxKAlphaDir= kAlphaDirs.get(j);
                    int k=Integer.valueOf( DxAxKAlphaDirName.split("_")[0].substring(1) );
                    float alpha=Float.valueOf( DxAxKAlphaDirName.split("_")[1].substring(1) );
                    System.out.println("db_build for k="+k+" alpha="+alpha+" in "+experimentLabel);
                    
                    StringBuilder sb=new StringBuilder();
                    sb.append("/usr/bin/java -jar "+exp.RAPPAJar.getAbsolutePath()+" ");
                    sb.append("-m b -a "+new Float(alpha).toString()+" -k "+new Integer(k).toString()+" ");
                    sb.append("-t "+Tx.getAbsolutePath()+" ");
                    sb.append("-i "+Ax.getAbsolutePath()+" ");
                    sb.append("-w "+DxAxKAlphaDir+" ");
                    sb.append("--ardir "+DxAxKAlphaDir+File.separator+"AR ");
                    sb.append("--extree "+DxAxKAlphaDir+File.separator+"extended_trees ");
                    sb.append("-v 1");

                    //execute script through qsub
                    sbQsubDBBUILDCommands.append(  "echo \""+sb.toString()+"\" |" +
                                            " qsub -N RAPx_"+experimentLabel+"_k"+k+"_a"+alpha+
                                            " -wd "+DxAxKAlphaDir.getAbsolutePath()+"\n");
                }
                
                
            }
            fw=new FileWriter(qsub_dbbuild_commands);
            fw.append(sbQsubDBBUILDCommands);
            fw.close();
            Files.setPosixFilePermissions(qsub_dbbuild_commands.toPath(), exp.perms);

            

            
            
//                   
            
            System.exit(1);
            


            
        } catch (IOException ex) {
            Logger.getLogger(RAPPASExperiment.class.getName()).log(Level.SEVERE, null, ex);
        } 

    }  
}
