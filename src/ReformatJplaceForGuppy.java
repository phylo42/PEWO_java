
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
            args=new String[4];
            args[0]="/home/benclaff/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/10000_aligned.jplace";
            args[1]="/home/benclaff/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/RAxML_portableTree.10000_aligned.jplace";
            args[2]="/home/benclaff/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/epa_result.jplace";
            args[3]="/home/benclaff/Dropbox/viromeplacer/test_datasets/KR_dist_experiment/placements_10000_last.fasta_union.jplace";
           
            
            System.out.println("ARGS: pplacerFile epaFile epangFile rappasFile");
            File pplacerFile=new File(args[0]);
            File epaFile=new File(args[1]);
            File epangFile=new File(args[2]);
            File rappasFile=new File(args[3]);
            
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
            fwJSON = new FileWriter(rappasFile.getParentFile().getAbsolutePath()+File.separator+"reform_RAPPAS.jplace");
            fwJSON.append(out);
            fwJSON.close();
            
            //------------------------------------------------
            //load epang file, update the tree, column order
            parser=new JSONParser();
            topLevel=(JSONObject)parser.parse(new FileReader(epangFile));
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
                                                         //old           new
                    valsNew.add(onePLacement.get(3));    //edge_num      distal_length  
                    valsNew.add(onePLacement.get(0));    //likelihood    edge_num  
                    valsNew.add(onePLacement.get(2));    //like_weight_r like_weight_r  
                    valsNew.add(onePLacement.get(1));    //distal_length likelihood  
                    valsNew.add(onePLacement.get(4));    //pendant_len   pendant_len
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
            System.out.println("EPANG: placementCount"+seqCount);
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
            fwJSON = new FileWriter(epangFile.getParentFile().getAbsolutePath()+File.separator+"reform_EPANG.jplace");
            fwJSON.append(out);
            fwJSON.close();
            
            
            
            //--------------------------------------------------------------
            //load epa file, update the tree, column order (multiplicity?)
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
                                                         //old           new
                    double d=(double)onePLacement.get(3);
                    if (d<0.001) {d=0;}
                    valsNew.add(d);                      //edge_num      distal_length  
                    valsNew.add(onePLacement.get(0));    //likelihood    edge_num  
                    valsNew.add(onePLacement.get(2));    //like_weight_r like_weight_r  
                    valsNew.add(onePLacement.get(1));    //distal_length likelihood  
                    valsNew.add(onePLacement.get(4));    //pendant_len   pendant_len
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
            System.out.println("EPA: placementCount"+seqCount);
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
            fwJSON = new FileWriter(epaFile.getParentFile().getAbsolutePath()+File.separator+"reform_EPA.jplace");
            fwJSON.append(out);
            fwJSON.close();
            
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(ReformatJplaceForGuppy.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ParseException ex) {
            Logger.getLogger(ReformatJplaceForGuppy.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    
}
