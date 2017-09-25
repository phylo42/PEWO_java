
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.PosixFilePermission;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Iterator;
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
public class RAPPASExperimentDBInRAM {
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
            String firstExperimentLabel=exp.prunedAlignmentsFiles.get(0).getName().split("\\.align$")[0];
            File DxDir=new File(exp.workDir.getCanonicalPath()+File.separator+"Dx");
            File DxA0Dir=new File(exp.workDir.getCanonicalPath()+File.separator+"Dx"+File.separator+firstExperimentLabel);
            File[] listFiles = DxA0Dir.listFiles();
            HashMap<Integer,Boolean> kmap=new HashMap<>();
            HashMap<Float,Boolean> alphamap=new HashMap<>();
            for (int i = 0; i < listFiles.length; i++) {
                File kAlphaDir = listFiles[i];
                if (kAlphaDir.getName().startsWith("k")) {
                    int k=Integer.parseInt(kAlphaDir.getName().split("_")[0].substring(1));
                    float alpha=Float.parseFloat(kAlphaDir.getName().split("_")[1].substring(1));
                    kmap.put(k, true);
                    alphamap.put(alpha, true);
                }
            }
            //convert to list to order 
            List<Integer> allK=kmap.keySet().stream().sorted().collect(Collectors.toList());
            List<Float> allAlpha=alphamap.keySet().stream().sorted().collect(Collectors.toList());
            
            //make a dir for qsub logs in Dx
            File qSubLogs=new File(DxDir.getAbsolutePath()+File.separator+"qsub_logs");
            qSubLogs.mkdir();
            

            
            ////////////////////////////////////////////////////////////////
            //build placement commands
            
            //file and filewriters holding list of commands
            File fileCommandList=null;
            FileWriter fwCommandList=null;
            

            
            //intermediate array loader scripts
            //version 18G (k<10)
            File shScript=new File(DxDir.getAbsolutePath()+File.separator+"array_loader_placement_18G.sh");
            //version 26G (k>=10)
            File shScript2=new File(DxDir.getAbsolutePath()+File.separator+"array_loader_placement_26G.sh");

            //then write the qsub_command in a file
            //ex: qsub -t 1-100 ./qsub_array_loader.sh /ngs/linard/tests_accuracy/pplacer_16s/Dx 
            File qsubCommand=new File(DxDir.getAbsolutePath()+File.separator+"qsub_rappas_placement");
            FileWriter fwQsubCommand=new FileWriter(qsubCommand);
            
            
            //For each k
            for (Iterator<Integer> iterator = allK.iterator(); iterator.hasNext();) {
                Integer k = iterator.next();
                
                //placement_commands files opened
                fileCommandList=new File(DxDir.getAbsolutePath()+File.separator+"k"+k+"_placement_commands.list");
                fwCommandList=new FileWriter(fileCommandList);

                //number of commands per k qsub array
                int commandCount=0;
                
                //list of rappas placement commands for this k
                StringBuilder sbPLACEMENTCommands=new StringBuilder();
            
                //for all AxDir, each combination of k/alpha
                for (int i = 0; i < exp.prunedAlignmentsFiles.size(); i++) {
                    File Ax = exp.prunedAlignmentsFiles.get(i);
                    File Tx = exp.prunedTreesFiles.get(i);
                    String experimentLabel=Ax.getName().split("\\.align$")[0];
                    //list Rx files corresponding to this Ax
                    File RxDir=new File(exp.workDir.getCanonicalPath()+File.separator+"Rx");
                    listFiles = RxDir.listFiles((File dir, String name) -> name.startsWith("R"+experimentLabel.substring(1))); //(A->R)XX_nxXXX
                    System.out.println("list"+Arrays.toString(listFiles));
                    
                    //for all alpha
                    for (Iterator<Float> iterator1 = allAlpha.iterator(); iterator1.hasNext();) {
                        Float alpha = iterator1.next();

                        String  kAlphaDirName= "k"+k+"_a"+alpha;
                        File  DxAxKAlphaDir= new File(DxDir.getAbsolutePath()+File.separator+experimentLabel+File.separator+kAlphaDirName);
                        System.out.println("placement for k="+k+" alpha="+alpha+" in "+experimentLabel+"/"+DxAxKAlphaDir.getName());

                        //for all query reads
                        for (int j = 0; j < listFiles.length; j++) {
                            File RxFile = listFiles[j];
                            //single, place on DB medium, then DB small
                            StringBuilder sb=new StringBuilder();
                            sb.append("/usr/bin/java ");
                            if (k<10)
                                sb.append("-Xmx16000m -Xms16000m ");
                            else {
                                sb.append("-Xmx24000m -Xms24000m ");
                            }
                            sb.append("-jar "+exp.RAPPAJar.getAbsolutePath()+" ");
                            sb.append("-XX:+UseG1GC -XX:ConcGCThreads=1 -XX:+PrintCommandLineFlags -XX:+AggressiveOpts "); //memory fine control
                            sb.append("-m b -a "+new Float(alpha).toString()+" -k "+new Integer(k).toString()+" ");
                            sb.append("-t "+Tx.getAbsolutePath()+" ");
                            sb.append("-i "+Ax.getAbsolutePath()+" ");
                            sb.append("-w "+DxAxKAlphaDir+" ");
                            sb.append("--ardir "+DxAxKAlphaDir+File.separator+"AR ");
                            //sb.append("--extree "+DxAxKAlphaDir+File.separator+"extended_trees "); //TODO: currently raise bugs, ids mappings problems...
                            sb.append("-v 1 ");
                            sb.append("--skipdbfull ");
                            sb.append("--froot "); //not necessary, as trees from Tx directory should already be rooted and with added_root node
                            sb.append("--dbinram ");
                            sb.append("-q "+RxFile.getAbsolutePath()+" ");
                            sb.append("--nsbound -100000.0 "); //skip calibration for this test.
                            sb.append("--unihash");
                            sb.append("\n");
                            sbPLACEMENTCommands.append(sb.toString());  

                            commandCount++;
                        }
                    }
                } 
                
                
                //placement_commands files closed
                fwCommandList.append(sbPLACEMENTCommands);
                fwCommandList.close();
                Files.setPosixFilePermissions(fileCommandList.toPath(), exp.perms);

                //add line in qsub_command file
                fwQsubCommand.append("qsub -pe smp 2 -N RAP_k"+k+" -t 1-"+commandCount+" -e "+qSubLogs.getAbsolutePath()+" -o "+qSubLogs.getAbsolutePath());
                if (k<8) {
                    fwQsubCommand.append(" "+shScript.getAbsolutePath());
                } else {
                    fwQsubCommand.append(" "+shScript2.getAbsolutePath());
                }
                fwQsubCommand.append(" "+DxDir.getAbsolutePath()+" "+fileCommandList.getName()+"\n");
                
                

            }
        
            //close qsub_command file
            fwQsubCommand.close();            
            Files.setPosixFilePermissions(qsubCommand.toPath(), exp.perms);
            
            //then write qsub array shell script in the working directory
            InputStream resourceAsStream = exp.getClass().getClassLoader().getResourceAsStream("scripts/array_loader_placement_18G.sh");
            Files.copy(resourceAsStream, shScript.toPath(), StandardCopyOption.REPLACE_EXISTING);
            resourceAsStream.close();
            Files.setPosixFilePermissions(shScript.toPath(), exp.perms);

            resourceAsStream = exp.getClass().getClassLoader().getResourceAsStream("scripts/array_loader_placement_26G.sh");
            Files.copy(resourceAsStream, shScript2.toPath(), StandardCopyOption.REPLACE_EXISTING);
            resourceAsStream.close();
            Files.setPosixFilePermissions(shScript2.toPath(), exp.perms);            
            
            //write script helping to relaunch failed jobs
            File shScript3=new File(DxDir.getAbsolutePath()+File.separator+"script_relaunch_failed.sh");
            resourceAsStream = exp.getClass().getClassLoader().getResourceAsStream("scripts/script_relaunch_failed.sh");
            Files.copy(resourceAsStream, shScript3.toPath(), StandardCopyOption.REPLACE_EXISTING);
            resourceAsStream.close();
            Files.setPosixFilePermissions(shScript3.toPath(), exp.perms);   

            System.out.println("DONE!");
            
            System.exit(0);
            


            
        } catch (IOException ex) {
            Logger.getLogger(RAPPASExperiment.class.getName()).log(Level.SEVERE, null, ex);
        } 

    }  
}
