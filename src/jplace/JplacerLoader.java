package jplace;


import etc.Infos;
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
    String treeString=null;
    
    //map containing the best placement
    //these are translated from jplace {x} nodeIds to the PhyloTree nodeIds
    //by this loader, i.e. this map contains the nodeIds of the PhyloTree,
    //not the ids of the jplace file itself
    //map(sequence_label)=ArrayList<JSONArray>
    HashMap<String,ArrayList<Integer>> nodeIds=new HashMap<>();
    HashMap<String,ArrayList<Double>> weightRatio=new HashMap<>();
    
    
    int edgeIdIndex=-1;
    int weightRatioIndex=-1;
    
    /**
     *
     * @param jplaceFile the value of jplaceFile
     * @param correctEPANGUnrooting the value of correctEPANGUnrooting
     */
    public JplacerLoader(File jplaceFile, boolean correctEPANGUnrooting) {
    
        try {
            
            
            JSONParser parser=new JSONParser();
            JSONObject topLevel=(JSONObject)parser.parse(new FileReader(jplaceFile));
            
            //read tree
            this.treeString=(String)topLevel.get("tree");
            
            //if this is from a EPANG output, a rooted input will be unrooted 
            // ((A,B),C)root; --> (C,B,A);
            //so, need to reorder string elements, then reroot the tree
            //as jplaceEdgeIds are registered in PhyloNodes, the rooted tree copy
            //returned by this operation still contains the correct jplace edge ids.
            //indexation is done after rooting.
            //jplaceId of added_root will be -1
            
            if (correctEPANGUnrooting) {
                Infos.println("Reversing EPA-ng unrooting...");

                //first change newick string to get root sons 
                //order from (C3,C2,C1); to (C1,C2,C3);

                //1st, define root
                int cladeClosingIndex=-1;
                for (int i = treeString.length()-1; i >= 0; i--) {
                    if (treeString.charAt(i)==')') {
                        cladeClosingIndex=i; 
                        break;
                    }
                }
                //System.out.println("Closing index: "+cladeClosingIndex);

                //2nd, extract C1 to C3
                String[] clades=new String[4]; //the 4th contains the root node 
                int depth=0;
                int cladeStart=1;
                int cladeCounter=0;
                for (int i = 0; i < treeString.length(); i++) {
                    char c=treeString.charAt(i);
                    if (c=='(') { depth++; }
                    if (c==')') { depth--; }
                    //we arrive or return to clade attached to root
                    if ( (depth==1 && c==',') || (depth==0 && i==cladeClosingIndex) ) {
                        //store previous clade string
                        if (i>0) {
                            clades[cladeCounter]=treeString.substring(cladeStart,i);
                            //System.out.println(clades[cladeCounter]);
                            cladeCounter++;
                        }
                        cladeStart=i+1;

                        continue;
                    }

                }
                //last one = the root
                clades[cladeCounter]=treeString.substring(cladeStart,treeString.length());
                //System.out.println(clades[cladeCounter]); 

                //reorder
                treeString="("+clades[2]+","+clades[1]+","+clades[0]+")"+clades[3];
                //System.out.println("---\n"+originalTreeString+"\n---\n");


                //then force rooting, which goes back to the original extended tree rooting
                this.tree=NewickReader.parseNewickTree2(treeString, true, true);

            } else {
                this.tree=NewickReader.parseNewickTree2(treeString, false, true); //init indexes already done in the parser
            }
            
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
                    int nodeId=tree.getJplaceMappingJPToNodeID(edgeJPlaceId.intValue());
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
    
    /**
     * return the tree as a newick string
     * @return 
     */
    public String getTreeString() {
        return treeString;
    }
    
    
    
    public static void main(String[] args) {
        //tests
        File test=new File("/home/ben/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML/logs/placements_mod_p4z1r36_query_only2.fasta_union.jplace");
        JplacerLoader jl=new JplacerLoader(test, false);

    }
    
    
    
    
}
