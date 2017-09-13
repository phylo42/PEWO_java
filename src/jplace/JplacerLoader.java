package jplace;


import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Stream;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.JSONValue;
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
    //map(sequence_label)=node_Id
    HashMap<String,Integer> bestPlacements=new HashMap<>();
    

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
            int edgeIdIndex=-1;
            JSONArray fields=(JSONArray)topLevel.get("fields");
            for (int i = 0; i < fields.size(); i++) {
                String field = (String)fields.get(i);
                if (field.equals("edge_num")) {
                    edgeIdIndex=i;
                    break;
                }
            }
            
            //load all placements
            JSONArray placements=(JSONArray)topLevel.get("placements");
            
            for (int i = 0; i < placements.size(); i++) {
                JSONObject pObject= (JSONObject)placements.get(i);
                //take best placement (1st in list)
                JSONArray pFields=(JSONArray)((JSONArray)pObject.get("p")).get(0);
                //System.out.println(pFields);
                //best edge
                //System.out.println(pFields.get(edgeIdIndex)+" "+pFields.get(edgeIdIndex).getClass());
                Long bestEdgeJPlaceId=(Long)pFields.get(edgeIdIndex);
                
                //list of sequences associate to this placement
                //this will be either a "n" field (name only,EPA) or a "nm" field (name mulitplicity,  pplacer and us)
                boolean isNm=true;
                JSONArray pNm=(JSONArray)pObject.get("nm");
                if (pNm==null) {
                    isNm=false;
                    pNm=(JSONArray)pObject.get("n");
                }
                //System.out.println("isNm:"+isNm);

                for (int j = 0; j < pNm.size(); j++) {
                    String n=null;
                    if (isNm) {     
                        n=(String)((JSONArray)pNm.get(j)).get(0);
                    } else {
                        n=(String)pNm.get(j);
                    }
                    //equivalent in phylotree
                    //System.out.println("bestEdgeJPlaceId:"+bestEdgeJPlaceId);
                    int nodeId=tree.getJplaceMapping(bestEdgeJPlaceId.intValue());
                    bestPlacements.put(n, nodeId);
                    //System.out.println( n+"->jplace:"+bestEdgeJPlaceId+"->phylotree_nodeId:"+nodeId
                    //                    +"->phylonode:"+tree.getById(nodeId));
                    
                }
                
                
            }

            
        } catch (IOException ex) {
            Logger.getLogger(JplacerLoader.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ParseException ex) {
            Logger.getLogger(JplacerLoader.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    
    }

    /**
     * map(sequence_name)=this.phylotree_nodeId
     * @return 
     */
    public HashMap<String, Integer> getBestPlacements() {
        return bestPlacements;
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
        File test=new File("/home/ben/Downloads/placementsR0_nx110_la_r300.fasta.jplace");
        File test2=new File("/home/ben/Downloads/placementsR0_nx110_la_r150.fasta.jplace");
        JplacerLoader jl=new JplacerLoader(test);
        System.out.println("##############");
        JplacerLoader jl2=new JplacerLoader(test2);
    }
    
    
    
    
}
