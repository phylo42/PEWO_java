package jplace;


import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
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
public class JplacerLoader {
    
    PhyloTree tree=null;
    
    //map containing the best placement
    //these are translated from jplace {x} nodeIds to the PhyloTree nodeIds
    //by this loader, i.e. this map contains the nodeIds of the PhyloTree,
    //not the ids of the jplace file itself
    //map(sequence_label)=ArrayList<JSONArray>
    HashMap<String,ArrayList<Integer>> nodeIds=new HashMap<>();
    HashMap<String,ArrayList<Double>> weightRatio=new HashMap<>();
    
    
    int edgeIdIndex=-1;
    int weightRatioIndex=-1;
    
    public JplacerLoader(File jplaceFile) {
    
        try {
            
            
            JSONParser parser=new JSONParser();
            JSONObject topLevel=(JSONObject)parser.parse(new FileReader(jplaceFile));
            
            //read tree
            String treeString=(String)topLevel.get("tree");
            tree=NewickReader.parseNewickTree2(treeString, false, true); //init indexes already done in the parser
            //System.out.println("isFromJplace:"+tree.isFromJplace());
            //System.out.println("jPlaceMappings:"+tree.getAllJPlaceMappings());
            
            
            //determine which p object columns are the edgeIds
            JSONArray fields=(JSONArray)topLevel.get("fields");
            for (int i = 0; i < fields.size(); i++) {
                String field = (String)fields.get(i);
                if (field.equals("edge_num")) {
                    edgeIdIndex=i;
                }
                if (field.equals("like_weight_ratio")) {
                    weightRatioIndex=i;
                }
            }
            
            //load all nodeIds
            JSONArray placements=(JSONArray)topLevel.get("placements");
            
            for (int i = 0; i < placements.size(); i++) {
                JSONObject placementsObject= (JSONObject)placements.get(i);
                //take best placement (1st in list)
                JSONArray pFields=(JSONArray)placementsObject.get("p");
                //System.out.println(pFields);
                //best edge
                //System.out.println(pFields.get(edgeIdIndex)+" "+pFields.get(edgeIdIndex).getClass());
                
                ArrayList<Integer> nodeIdList=new ArrayList<>(pFields.size());
                ArrayList<Double> weightRatiosList=new ArrayList<>(pFields.size());
                for (int j = 0; j < pFields.size(); j++) {
                    JSONArray stats = (JSONArray)pFields.get(j);
                    Long edgeJPlaceId=(Long)stats.get(edgeIdIndex);
                    //equivalent in phylotree
                    int nodeId=tree.getJplaceMapping(edgeJPlaceId.intValue());
                    nodeIdList.add(nodeId);                    
                    Object o=stats.get(weightRatioIndex);
                    try {
                        if (o.getClass().getName().contains("Double")) {
                            weightRatiosList.add((Double)stats.get(weightRatioIndex));
                        } else if (o.getClass().getName().contains("Long")) {
                            weightRatiosList.add(new Double(((Long)stats.get(weightRatioIndex)).doubleValue()));
                        }
                    } catch (Exception ex) {
                        ex.printStackTrace();
                        if (o==null) {
                            System.out.println("Weight Ratio parsed as 'null'.");
                        }
                        System.out.println("Line: "+stats.toString());
                        weightRatiosList.add(0.0);
                        
                    }
                    
                   
                }
                
                
                
                //list of sequences associate to this placement
                //this will be either a "n" field (name only,EPA) or a "nm" field (name mulitplicity,  pplacer and us)
                boolean isNm=true;
                JSONArray pNm=(JSONArray)placementsObject.get("nm");
                if (pNm==null) {
                    isNm=false;
                    pNm=(JSONArray)placementsObject.get("n");
                }
                //System.out.println("isNm:"+isNm);

                for (int j = 0; j < pNm.size(); j++) {
                    String n=null;
                    if (isNm) {     
                        n=(String)((JSONArray)pNm.get(j)).get(0);
                    } else {
                        n=(String)pNm.get(j);
                    }
                    this.nodeIds.put(n, nodeIdList);
                    this.weightRatio.put(n, weightRatiosList);
                }
                
                
            }

            
        } catch (IOException ex) {
            Logger.getLogger(JplacerLoader.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ParseException ex) {
            Logger.getLogger(JplacerLoader.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    
    }

    /**
     * map(sequence_name)=this.phylotree_nodeIds
     * @return 
     */
    public HashMap<String, ArrayList<Integer>> getNodeIds() {
        return nodeIds;
    }
    
    /**
     * map(sequence_name)=weightRatios
     * @return 
     */
    public HashMap<String, ArrayList<Double>> getWeightRatios() {
        return weightRatio;
    }

    /**
     * return phylotree equivalent of the newick stored in the jplace
     * @return 
     */
    public PhyloTree getTree() {
        return tree;
    }
    
    
    
    
    
    public static void main(String[] args) {
        //tests
        File test=new File("/home/ben/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML/logs/placements_mod_p4z1r36_query_only2.fasta_union.jplace");
        JplacerLoader jl=new JplacerLoader(test);

    }
    
    
    
    
}
