/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dtx;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ben
 */
public class D2tx {

    //access to correct Dtx line through node Label
    HashMap<String, Integer> pruningLabels=new HashMap<>(); //map(prunedNodeLabel)=dataLineIndex
    //access to correct Dtx line through node id
    HashMap<Integer, Integer> pruningIds=new HashMap<>();  //map(prunedNodeId)=dataLineIndex

    HashMap<String,Integer> colLabels=new HashMap<>(); //map(label)=dataColIndex
    HashMap<Integer,Integer>  colIds=new HashMap<>();  //map(nodeId)=dataColIndex
    ArrayList<ArrayList<Double>> data=new ArrayList<>(); //data themselves


    public D2tx(File Dtx) {
        BufferedReader br=null;
        try {
            br = new BufferedReader(new FileReader(Dtx));
            String line=null;
            int lineIndex=0;
            data=new ArrayList<>();
            while ((line=br.readLine())!=null) {                
                String[] elts=line.split(";");
                //1st line
                if (elts[0].equals("nodeLabels")) {
                    for (int i = 2; i < elts.length; i++) {
                        colLabels.put(elts[i],i-2);
                    }
                    continue;
                }
                //2nd line
                if (elts[1].equals("nodeIds")) {
                    for (int i = 2; i < elts.length; i++) {
                        colIds.put(Integer.parseInt(elts[i]),i-2);
                    }
                    continue;
                }
                //other lines
                ArrayList<Double> lineData=new ArrayList<>();
                pruningLabels.put(elts[0], lineIndex);
                pruningIds.put(Integer.parseInt(elts[1]), lineIndex);
                for (int i = 2; i < elts.length; i++) {
                    lineData.add(Double.parseDouble(elts[i]));
                }
                data.add(lineData);
                lineIndex++;
            }
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(D2tx.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(D2tx.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(D2tx.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    /**
     * give the pruning experiment through Nx id, retrieve node distance
     * of nodeId to the edge where Nx leaves are expected to be placed
     * @param NxId
     * @param nodeId
     * @return -1 if nodeId is not anymore in this Nx pruning, the node distance otherwise
     */
    public Double getNodeDistance(int NxId,int nodeId) {
        int dataLine=pruningIds.get(NxId);
        int dataColumn=colIds.get(nodeId);
        return data.get(dataLine).get(dataColumn);
    }

    @Override
    public String toString() {
        StringBuilder sb=new StringBuilder();
        for (int i = 0; i < data.size(); i++) {
            ArrayList<Double> get = data.get(i);
            for (int j = 0; j < get.size(); j++) {
                Double get1 = get.get(j);
                if (j>0)
                    sb.append(";");
                sb.append(get1);
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
    
    
    
}
