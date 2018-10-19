
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.BufferedWriter;
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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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
public class EPANGExperiment {
    
    //workDir
    String HOME = System.getenv("HOME");
    File workDir=new File(HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/6_leaves_test_set");
    //set if analysis is protein or DNA/RNA
    boolean proteinAnalysis=false;

    //list of new files
    public List<File> prunedAlignmentsFiles=new ArrayList<>(); //list of pruned alignments
    public List<File> prunedTreesFiles=new ArrayList<>(); //list of pruned Trees
   
    //where is EPA 
    File EPANGBinary=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/epa_sse3/bin/epa-ng");

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
            System.out.println("ARGS: workDir EPANGBinary [nucl=0|prot=1]");
            
            //launch
            EPANGExperiment exp=new EPANGExperiment();
            //LOAD ALL EXPERIMENTS FOUND IN WORK DIR
            ///////////////////////////////////////////////////
            if(args.length>0) {
                exp.workDir=new File(args[0]);
                exp.EPANGBinary=new File(args[1]);
                int protein=Integer.parseInt(args[2]);
                exp.proteinAnalysis=(protein>0);
                System.out.println("workDir: "+exp.workDir);
                System.out.println("EPANGBinary: "+exp.EPANGBinary);
                System.out.println("proteinAnalysis:"+exp.proteinAnalysis);
            }  
            
            //pattern to distinguish queriesBuf from refsBuf in hmmalign results
            String pattern="_r[0-9]+_[0-9]+_[0-9]+(_[0-9]+)?$";
            Pattern p=Pattern.compile(pattern);

            //load alignments trees from Ax/Tx directories
            File Ax=new File(exp.workDir.getAbsolutePath()+File.separator+"Ax");
            File Tx=new File(exp.workDir.getAbsolutePath()+File.separator+"Tx");
            exp.prunedAlignmentsFiles=Arrays.stream(Ax.listFiles()).sorted().collect(Collectors.toList());
            exp.prunedTreesFiles=Arrays.stream(Tx.listFiles()).filter((f)->f.getName().startsWith("T")).sorted().collect(Collectors.toList());
            System.out.println(exp.prunedAlignmentsFiles);
            System.out.println(exp.prunedTreesFiles);
            if (    exp.prunedAlignmentsFiles.size()<1 ||
                    exp.prunedTreesFiles.size()<1 ||
                    exp.prunedAlignmentsFiles.size()!=exp.prunedTreesFiles.size()) {
                System.out.println("Unexpected Ax/Tx files");
            }   //build directory /RAXMLx
            File EPAxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"EPANGx");
            EPAxDir.mkdir();
            File HMMxDir=new File(exp.workDir.getAbsolutePath()+File.separator+"HMMx");
            //list of raxml commands
            StringBuilder sbQsubCommands=new StringBuilder();
            //build all Ax directories in /HMMx and for each set of hmm-aligned
            //reads (/HMMx/Ax_nx_la/Rx_nx_la_r), build the list of commands
            for (int i=0;i<exp.prunedAlignmentsFiles.size();i++) {
                String experimentLabel=exp.prunedAlignmentsFiles.get(i).getName().split("\\.align$")[0];
                System.out.println("EPANG commands for "+experimentLabel);
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
                
                //create epa_pipeline.sh epaScript
                File epaScript=new File(EPAxAxDir.getAbsolutePath()+File.separator+"epang_pipeline.sh");
                BufferedWriter fwEpaScript = Files.newBufferedWriter(epaScript.toPath());
                //add command for each read
                for (int j = 0; j < listFiles.length; j++) {
                    File f = listFiles[j];
                    String readAlignLabel=f.getName().split("\\.aln\\.fasta$")[0];
                    
                    //epa-ng require to split reference and queriesBuf from the same 
                    //multiple alignment (need same number of sites)
                    //we load this fasta, and write 2 new fastas:
                    //one with queriesBuf, wich finishes with _r[0-9]+_[0-9]+_[0-9]+$
                    //one with refsBuf, which are the remaining ones
                    StringBuilder queriesBuf=new StringBuilder();
                    StringBuilder refsBuf=new StringBuilder();
                    FASTAPointer fp=new FASTAPointer(f, false);
                    Fasta seq=null;
                    while ((seq=fp.nextSequenceAsFastaObject())!=null) {
                        //regexp requires the end addition due to read duplicates (simulated reads experiment)
                        Matcher matcher = p.matcher(seq.getHeader());
                        if (matcher.find()) {
                            queriesBuf.append(seq.getFormatedFasta()+"\n");
                        } else {
                            refsBuf.append(seq.getFormatedFasta()+"\n");
                        }
                    }
                    fp.closePointer();
                    File queriesFile=new File(f.getAbsolutePath()+"_queries");
                    BufferedWriter queriesBW = Files.newBufferedWriter(queriesFile.toPath());
                    queriesBW.append(queriesBuf);
                    queriesBW.close();
                    File refsFile=new File(f.getAbsolutePath()+"_refs");
                    BufferedWriter refsBW = Files.newBufferedWriter(refsFile.toPath());
                    refsBW.append(refsBuf);
                    refsBW.close();
                    //build the placement commands
                    //example:  ../epa_sse3/bin/epa-ng -T 1 -t tree -s align -q queriesBuf --verbose ;
                    File statFile=new File(Tx.getAbsoluteFile()+File.separator+"RAxML_info.OPTIM_"+exp.prunedTreesFiles.get(i).getName());
                    StringBuilder sbRAxMLCommand=new StringBuilder();
                    //for parallel version
                    //sbRAxMLCommand.append(exp.EPANGBinary.getAbsolutePath()+" --verbose -T 1 ");
                    //EPANG was recompiled with make EPA_SERIAL=1, which removes MPI, then option -T is refused
                    sbRAxMLCommand.append(exp.EPANGBinary.getAbsolutePath()+" --verbose ");

                    sbRAxMLCommand.append(
                            "-w "+EPAxAxDir.getAbsolutePath()+" " +
                            "-q "+queriesFile.getAbsolutePath()+" " +
                            "-s "+refsFile.getAbsolutePath()+" " +
                            "-t "+exp.prunedTreesFiles.get(i).getAbsolutePath()+" "+
                            "-m "+statFile.getAbsolutePath()+" "+
                            "\n"
                    );
                    //current epa alpha version cannot rename outputs, need to be done manually
                    sbRAxMLCommand.append("mv "+
                                            EPAxAxDir.getAbsolutePath()+File.separator+"epa_info.log "+
                                            EPAxAxDir.getAbsolutePath()+File.separator+queriesFile.getName()+"_epa_info.log \n"
                                        );
                    sbRAxMLCommand.append("mv "+
                                            EPAxAxDir.getAbsolutePath()+File.separator+"epa_result.jplace "+
                                            EPAxAxDir.getAbsolutePath()+File.separator+queriesFile.getName()+"_epa_result.jplace \n"
                                        );
                    
                    fwEpaScript.append(sbRAxMLCommand.toString());
                    sbRAxMLCommand=null;
                    fwEpaScript.append("\n");
                }
                fwEpaScript.close();
                Files.setPosixFilePermissions(epaScript.toPath(), exp.perms);
                
                //add qsub epa script call 
                sbQsubCommands.append("echo \""+epaScript.getAbsolutePath()+"\" |  qsub -N EPANGx_"+experimentLabel+
                        " -wd "+EPAxAxDir.getAbsolutePath()+"\n");
            }
            File qsubEPACommands=new File(EPAxDir.getAbsolutePath()+File.separator+"qsub_epa_commands");
            fw = new FileWriter(qsubEPACommands);
            fw.append(sbQsubCommands.toString());
            fw.close();
            Files.setPosixFilePermissions(qsubEPACommands.toPath(), exp.perms);
            
        } catch (IOException ex) {
            Logger.getLogger(EPAExperiment.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(EPAExperiment.class.getName()).log(Level.SEVERE, null, ex);
            }  
        }
            

    }  
}
