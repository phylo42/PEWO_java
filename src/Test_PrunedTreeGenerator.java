
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class Test_PrunedTreeGenerator {
    
    public static void main(String[] args) {
        
        args=new String[16];
        args[0]="/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/greengenes";
        args[1]="/home/ben/Dropbox/viromeplacer/test_datasets/software/hmmer-3.1b2/binaries";
        args[2]="/home/ben/Dropbox/viromeplacer/test_datasets/software/RAxML-8.2.9/raxmlHPC-PTHREADS-SSE3";
        args[3]="/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/greengenes/align.reduced";
        args[4]="/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/greengenes/RAxML_result.OPTIM.REROOT_BACT-ARCH";
        args[5]="100";
        args[6]="150,300,600,1200";
        args[7]="15";
        args[8]="1";
        args[9]="6";
        args[10]="12";
        args[11]="2";
        args[12]="1.00";
        args[13]="1.00";
        args[14]="0.25";
        args[15]="0";
        
        PrunedTreeGenerator.main(args);
        
    }
  
}
