
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ListIterator;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import jplace.JplacerLoader;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import tree.PhyloNode;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author benclaff
 */
public class ReformatJplaceForGuppy {
    
    public static void main(String[] args) {
        
        FileWriter fwJSON=null;
        try {
            
            //DEBUG
            System.out.println("ARGS: PPL.jplace RAXML.jplace EPAPREPL.jplace EPATHOR.jplace RAPPASGUPPY.jplace");
//            args=new String[5];
//            args[0]="/home/ben/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/p4z1r36_aligned.jplace";
//            args[1]="/home/ben/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/RAxML_portableTree.p4z1r36_aligned.jplace";
//            args[2]="/home/ben/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/epa_result_preplacement.jplace";
//            args[3]="/home/ben/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/epang_result_thorough.jplace";
//            args[4]="/home/ben/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/placements_p4z1r36_query_only_no_gaps.fasta_union.jplace";
            
            
            File pplacerFile=new File(args[0]);
            File epaFile=new File(args[1]);
            File epangPreplacementFile=new File(args[2]);
            File epangThoroughFile=new File(args[3]);
            File rappasFile=new File(args[4]);
            
            NumberFormat nf=null;
            nf=NumberFormat.getNumberInstance(Locale.UK);
            nf.setMaximumFractionDigits(8);
            nf.setMinimumFractionDigits(8);
            nf.setParseIntegerOnly(false);

            //load pplacer filer, which will be the base
            //from which we extract the tree
            JplacerLoader jpl=new JplacerLoader(pplacerFile, false);
            String pplacerTree=jpl.getTreeString();
            
            //-----------------------------------------
            //load rappas file, only update the tree
            System.out.println("UPDATE RAPPAS -----------------------------------");
            JSONParser parser=new JSONParser();
            JSONObject topLevel=(JSONObject)parser.parse(new FileReader(rappasFile));
            //update tree
            topLevel.put("tree",pplacerTree);
            String out=topLevel.toJSONString();
            //just some basic formatting for human readability
            out=out.replaceAll("\\},\\{", "\n\\},\\{\n\t"); //},{
            out=out.replaceAll("\\],\"","\\],\n\t\"");   //],"
            out=out.replaceAll("\\]\\}\\],", "\\]\n\\}\n\\],\n"); //]}]
            out=out.replaceAll(",\"placements\":\\[\\{\"p\"", ",\n\"placements\":\n[\n{\n\t\"p\"");
            out=out.replaceAll("\\],\\[", "\\],\n\t\\[");
            out=out.replaceAll("\"p\":\\[\\[","\"p\":\n\t\\[\\[");
            out=out.replaceAll("\"nm\":\\[\\[","\"nm\":\n\t\\[\\[");
            //out=out.replace("]},", "]},"); //]}
            fwJSON = new FileWriter(rappasFile.getParentFile().getAbsolutePath()+File.separator+"re_RAPPAS.jplace");
            fwJSON.append(out);
            fwJSON.close();
            
            //------------------------------------------------
            //load epang file, update the tree, column order
            System.out.println("UPDATE EPANG-PREPLACEMENT -----------------------------------");
            parser=new JSONParser();
            topLevel=(JSONObject)parser.parse(new FileReader(epangPreplacementFile));
            //1.change metadata
            JSONArray fList=new JSONArray();
            fList.add("distal_length");
            fList.add("edge_num"); 
            fList.add("like_weight_ratio");
            fList.add("likelihood");
            fList.add("pendant_length");
            topLevel.put("fields", fList);
            //2. change all p objects: n to nm and add multiplicity of 1
            JSONArray placements = (JSONArray)topLevel.get("placements");
            int placementCount=0;
            int seqCount=0;
            for (ListIterator it=placements.listIterator();it.hasNext();) {
                placementCount++;
                JSONObject currentPlacement=(JSONObject)it.next();
                JSONArray pVals=(JSONArray)currentPlacement.get("p");
                JSONArray nVals=(JSONArray)currentPlacement.get("n");
                //iterate on value lists
                for (int i=0;i<pVals.size();i++) {
                    JSONArray onePLacement=(JSONArray)pVals.get(i);
                    JSONArray valsNew=new JSONArray();
                    
                    Double distal_length=(double)onePLacement.get(3);
                    //check is < to the branch length (issues related to 
                    //EPA which optimizes the tree whatever reference tree
                    //is given to him. This make that sometimes the 
                    //pendant_length is longer than the branch length
                    //which cause a guppy error
                    Long jplaceEdgeId=(long)onePLacement.get(0);
                    int nodeId=jpl.getTree().getAllJPlaceMappingsJPToNodeID().get(jplaceEdgeId.intValue());
                    PhyloNode node=jpl.getTree().getById(nodeId);
                    float bl=node.getBranchLengthToAncestor();
                    if (bl<distal_length.floatValue()) {
                        distal_length=bl-0.00000001;
                        if (distal_length<0) {distal_length=0.0;}
                        System.out.print("Test: distal_len= "+onePLacement.get(3)+"(label="+node.getLabel()+") VS "+jpl.getTree().getById(nodeId));
                        System.out.println(" ==> change distal_length="+distal_length);
                    }
                                                         //old           new
                    valsNew.add(distal_length);    //edge_num      distal_length  
                    valsNew.add(onePLacement.get(0));    //likelihood    edge_num  
                    valsNew.add(onePLacement.get(2));    //like_weight_r like_weight_r  
                    valsNew.add(onePLacement.get(1));    //distal_length likelihood  
                    valsNew.add(0.0);    //pendant_len   pendant_len
                    pVals.set(i, valsNew);
                }
                
                //remove "n", replace by "nm" with reordered columns
                currentPlacement.remove("n");
                JSONArray nmVals=new JSONArray();
                for (int i = 0; i < nVals.size(); i++) {
                    JSONArray id=new JSONArray();
                    id.add(nVals.get(i));
                    id.add(1);
                    nmVals.add(id);
                    seqCount++;
                }
                currentPlacement.put("nm", nmVals);
            }
            System.out.println("EPANG: seqCount"+seqCount);
            System.out.println("EPANG: placementCount"+placementCount);
            //3. change tree
            topLevel.put("tree", pplacerTree);
            //4. output
            out=topLevel.toJSONString();
            //out=out.replaceAll("[0-9\\.]+E-[0-9]+", "0.0001"); //HACK necessary because scientific numbers shows issues in Guppy
            out=out.replaceAll("\\},\\{", "\n\\},\\{\n\t"); //},{
            out=out.replaceAll("\\],\"","\\],\n\t\"");   //],"
            out=out.replaceAll("\\]\\}\\],", "\\]\n\\}\n\\],\n"); //]}]
            out=out.replaceAll(",\"placements\":\\[\\{\"p\"", ",\n\"placements\":\n[\n{\n\t\"p\"");
            out=out.replaceAll("\\],\\[", "\\],\n\t\\[");
            out=out.replaceAll("\"p\":\\[\\[","\"p\":\n\t\\[\\[");
            out=out.replaceAll("\"nm\":\\[\\[","\"nm\":\n\t\\[\\[");
            fwJSON = new FileWriter(epangPreplacementFile.getParentFile().getAbsolutePath()+File.separator+"re_EPANG-PREPLACEMENT.jplace");
            fwJSON.append(out);
            fwJSON.close();
            
            //------------------------------------------------
            //load epang file, update the tree, column order
            System.out.println("UPDATE EPANG-THOROUGH -----------------------------------");
            parser=new JSONParser();
            topLevel=(JSONObject)parser.parse(new FileReader(epangThoroughFile));
            //1.change metadata
            fList=new JSONArray();
            fList.add("distal_length");
            fList.add("edge_num"); 
            fList.add("like_weight_ratio");
            fList.add("likelihood");
            fList.add("pendant_length");
            topLevel.put("fields", fList);
            //2. change all p objects: n to nm and add multiplicity of 1
            placements = (JSONArray)topLevel.get("placements");
            placementCount=0;
            seqCount=0;
            for (ListIterator it=placements.listIterator();it.hasNext();) {
                placementCount++;
                JSONObject currentPlacement=(JSONObject)it.next();
                JSONArray pVals=(JSONArray)currentPlacement.get("p");
                JSONArray nVals=(JSONArray)currentPlacement.get("n");
                //iterate on value lists
                for (int i=0;i<pVals.size();i++) {
                    JSONArray onePLacement=(JSONArray)pVals.get(i);
                    JSONArray valsNew=new JSONArray();
                    
                    Double distal_length=(double)onePLacement.get(3);
                    //check is < to the branch length (issues related to 
                    //EPA which optimizes the tree whatever reference tree
                    //is given to him. This make that sometimes the 
                    //pendant_length is longer than the branch length
                    //which cause a guppy error
                    Long jplaceEdgeId=(long)onePLacement.get(0);
                    int nodeId=jpl.getTree().getAllJPlaceMappingsJPToNodeID().get(jplaceEdgeId.intValue());
                    PhyloNode node=jpl.getTree().getById(nodeId);
                    float bl=node.getBranchLengthToAncestor();
                    if (bl<distal_length.floatValue()) {
                        distal_length=bl-0.00000001;
                        if (distal_length<0) {distal_length=0.0;}
                        System.out.print("Test: distal_len= "+onePLacement.get(3)+"(label="+node.getLabel()+") VS "+jpl.getTree().getById(nodeId));
                        System.out.println(" ==> change distal_length="+distal_length);
                    }
                                                         //old           new
                    valsNew.add(distal_length);    //edge_num      distal_length  
                    valsNew.add(onePLacement.get(0));    //likelihood    edge_num  
                    valsNew.add(onePLacement.get(2));    //like_weight_r like_weight_r  
                    valsNew.add(onePLacement.get(1));    //distal_length likelihood  
                    valsNew.add(0.0);    //pendant_len   pendant_len
                    pVals.set(i, valsNew);
                }
                
                //remove "n", replace by "nm" with reordered columns
                currentPlacement.remove("n");
                JSONArray nmVals=new JSONArray();
                for (int i = 0; i < nVals.size(); i++) {
                    JSONArray id=new JSONArray();
                    id.add(nVals.get(i));
                    id.add(1);
                    nmVals.add(id);
                    seqCount++;
                }
                currentPlacement.put("nm", nmVals);
            }
            System.out.println("EPANG-THOROUGH: seqCount"+seqCount);
            System.out.println("EPANG-THOROUGH: placementCount"+placementCount);
            //3. change tree
            topLevel.put("tree", pplacerTree);
            //4. output
            out=topLevel.toJSONString();
            //out=out.replaceAll("[0-9\\.]+E-[0-9]+", "0.0001"); //HACK necessary because scientific numbers shows issues in Guppy
            out=out.replaceAll("\\},\\{", "\n\\},\\{\n\t"); //},{
            out=out.replaceAll("\\],\"","\\],\n\t\"");   //],"
            out=out.replaceAll("\\]\\}\\],", "\\]\n\\}\n\\],\n"); //]}]
            out=out.replaceAll(",\"placements\":\\[\\{\"p\"", ",\n\"placements\":\n[\n{\n\t\"p\"");
            out=out.replaceAll("\\],\\[", "\\],\n\t\\[");
            out=out.replaceAll("\"p\":\\[\\[","\"p\":\n\t\\[\\[");
            out=out.replaceAll("\"nm\":\\[\\[","\"nm\":\n\t\\[\\[");
            fwJSON = new FileWriter(epangPreplacementFile.getParentFile().getAbsolutePath()+File.separator+"re_EPANG-THOROUGH.jplace");
            fwJSON.append(out);
            fwJSON.close();
            
            
            
            //--------------------------------------------------------------
            //load epa file, update the tree, column order (multiplicity?)
            System.out.println("UPDATE EPA -----------------------------------");
            parser=new JSONParser();
            topLevel=(JSONObject)parser.parse(new FileReader(epaFile));
            //1.change metadata
            fList=new JSONArray();
            fList.add("distal_length");
            fList.add("edge_num"); 
            fList.add("like_weight_ratio");
            fList.add("likelihood");
            fList.add("pendant_length");
            topLevel.put("fields", fList);
            topLevel.put("version", 3);
            //2. change all p objects: n to nm and add multiplicity of 1
            seqCount=0;
            placementCount=0;
            placements = (JSONArray)topLevel.get("placements");
            for (ListIterator it=placements.listIterator();it.hasNext();) {
                placementCount++;
                JSONObject currentPlacement=(JSONObject)it.next();
                JSONArray pVals=(JSONArray)currentPlacement.get("p");
                JSONArray nVals=(JSONArray)currentPlacement.get("n");
                //iterate on value lists
                for (int i=0;i<pVals.size();i++) {
                    JSONArray onePLacement=(JSONArray)pVals.get(i);
                    JSONArray valsNew=new JSONArray();
                    
                    Double distal_length=(double)onePLacement.get(3);
                    //check is < to the branch length (issues related to 
                    //EPA which optimizes the tree whatever reference tree
                    //is given to him. This make that sometimes the 
                    //pendant_length is longer than the branch length
                    //which cause a guppy error
                    Long jplaceEdgeId=(long)onePLacement.get(0);
                    int nodeId=jpl.getTree().getAllJPlaceMappingsJPToNodeID().get(jplaceEdgeId.intValue());
                    float bl=jpl.getTree().getById(nodeId).getBranchLengthToAncestor();
                    if (bl<distal_length.floatValue()) {
                        distal_length=bl-0.000001;
                        if (distal_length<0) {distal_length=0.0;}
                        System.out.print("Test: distal_len= "+onePLacement.get(3)+" VS "+jpl.getTree().getById(nodeId));
                        System.out.println(" ==> change distal_length="+distal_length);
                    }
                                                         //old           new
                    valsNew.add(distal_length);          //edge_num      distal_length  
                    valsNew.add(onePLacement.get(0));    //likelihood    edge_num  
                    valsNew.add(onePLacement.get(2));    //like_weight_r like_weight_r  
                    valsNew.add(onePLacement.get(1));    //distal_length likelihood  
                    valsNew.add(0.0);    //pendant_len   pendant_len
                    pVals.set(i, valsNew);
                }
                //remove "n", replace by "nm" with reordered columns
                currentPlacement.remove("n");
                JSONArray nmVals=new JSONArray();
                for (int i = 0; i < nVals.size(); i++) {
                    JSONArray id=new JSONArray();
                    id.add(nVals.get(i));
                    id.add(1);
                    nmVals.add(id);
                    seqCount++;
                }
                currentPlacement.put("nm", nmVals);
            }
            System.out.println("EPA: seqCount"+seqCount);
            System.out.println("EPA: placementCount"+placementCount);
            //3. change tree
            topLevel.put("tree", pplacerTree);
            //4. output
            out=topLevel.toJSONString();
            //out=out.replaceAll("[0-9\\.]+E-[0-9]+", "0.0001"); //HACK necessary because scientific numbers shows issues in Guppy
            out=out.replaceAll("\\},\\{", "\n\\},\\{\n\t"); //},{
            out=out.replaceAll("\\],\"","\\],\n\t\"");   //],"
            out=out.replaceAll("\\]\\}\\],", "\\]\n\\}\n\\],\n"); //]}]
            out=out.replaceAll(",\"placements\":\\[\\{\"p\"", ",\n\"placements\":\n[\n{\n\t\"p\"");
            out=out.replaceAll("\\],\\[", "\\],\n\t\\[");
            out=out.replaceAll("\"p\":\\[\\[","\"p\":\n\t\\[\\[");
            out=out.replaceAll("\"nm\":\\[\\[","\"nm\":\n\t\\[\\[");
            fwJSON = new FileWriter(epaFile.getParentFile().getAbsolutePath()+File.separator+"re_EPA.jplace");
            fwJSON.append(out);
            fwJSON.close();
            
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(ReformatJplaceForGuppy.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ParseException ex) {
            Logger.getLogger(ReformatJplaceForGuppy.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    
}
