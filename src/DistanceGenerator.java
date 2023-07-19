
import dtx.Dtx;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import jplace.JplacerLoader;
import jplace.JplacerLoader.Placement;
import tree.PhyloTree;

/**
 * Compute the Node Distance (ND) and expected Node Distance (eND)
 * using all jplace files found in the PEWO_wordir/SOFTWARE_NAME/ directories
 * Note that:
 * 1) the actual node distance between all tree node pairs
 * are pre-computed at pruning step and saved as a matrix (DTx file)
 * 2) the expected placement is saved in expected_placements.bin
 *    (which contains not only this data but also the full pruned tree, trifurcations...)
 *
 * @author benjamin linard
 */
public class DistanceGenerator {
    
    File workDir;
    ArrayList<File> dirs;
    LinkedHashMap<String,Integer> paramSet = new LinkedHashMap<>(32);

    // Pattern: ([0-9]+)_r([0-9]+)_([a-z]+)([0-9A-Z]+)_.*_([a-z]+)\.jplace
    final Pattern patternMatch = Pattern.compile("([a-z]+)([0-9A-z\\-\\.]+)");

    DistanceGenerator(File workingDirectory, ArrayList<File> softwareDirs) {
        workDir = workingDirectory;
        dirs = softwareDirs;
    }

    /// Creates the set of parameters evaluated by the workflow
    ///     "parameter" -> "integer"
    LinkedHashMap<String,Integer> initializeParamSet(List<Path> jplaceFiles) {
        int columnCounter = initializeCommonColumns();
        initializeUncommonColumns(jplaceFiles, columnCounter);
        // System.out.println("Param set: "+paramSet);
        return paramSet;
    }

    /// Initializes the paramSet with the columns common to all software
    int initializeCommonColumns() {
        int columnCounter = -1;
        paramSet.put("software", ++columnCounter);
        paramSet.put("pruning", ++columnCounter);
        paramSet.put("query", ++columnCounter);
        paramSet.put("rstart", ++columnCounter);
        paramSet.put("rend", ++columnCounter);
        paramSet.put("nd", ++columnCounter);
        paramSet.put("e_nd", ++columnCounter);
        return columnCounter;
    }

    /// Updates the paramSet with the columns specific to tested
    /// software based on jplace file names
    void initializeUncommonColumns(List<Path> jplaceFiles, int columnCounter) {
        for (Path p:jplaceFiles) {
            String filename = p.toFile().getName();
            String[] globalGroups = filename.split("_");

            // 1st element is pruning id, last element is "software.jplace"
            for (int i = 1; i < globalGroups.length - 1; i++) {
                Matcher m = patternMatch.matcher(globalGroups[i]);
                if (m.matches()) {
                    String param = m.group(1);
                    if (!paramSet.containsKey(param)) {
                        paramSet.put(m.group(1), ++columnCounter);
                    }
                }
            }
        }
    }

    // Finds all .jplace files in the list of directories
    static List<Path> findJplaceFiles(List<File> placementSoftwareDirs) throws IOException {
        List<Path> allJplaceFiles = new ArrayList<>();
        for (File dir : placementSoftwareDirs) {
            System.out.println("Scanning for "+dir.getName()+" jplace results...");
            List<Path> jplaceFiles = Files.find(
                    dir.toPath(),
                    999,
                    (p,b)-> b.isRegularFile() && p.getFileName().toString().endsWith(".jplace")
            ).collect(Collectors.toList());
            System.out.println("# jplace found: "+jplaceFiles.size());
            allJplaceFiles.addAll(jplaceFiles);
        }

        return allJplaceFiles;
    }

    void computeDistances() throws IOException, ClassNotFoundException {

        // LOAD BINARY INFORMATION
        ////////////////////////////////////////////////////

        // load the expected placements
        File expPLaceFile = new File(workDir + File.separator + "expected_placements.bin");
        System.out.println("Loading "+expPLaceFile.getAbsolutePath());
        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(expPLaceFile),4096));
        System.out.println("Loading NxIndex");
        HashMap<Integer,Integer> NxIndex = (HashMap<Integer,Integer>)ois.readObject();
        System.out.println("Loading pruningIndex");
        HashMap<Integer,Integer> pruningIndex = (HashMap<Integer,Integer>)ois.readObject();
        System.out.println("Loading expected placements");
        ArrayList<ArrayList<Integer>> expectedPlacements = (ArrayList<ArrayList<Integer>>)ois.readObject();
        System.out.println("Loading trees");
        ArrayList<PhyloTree> experimentTrees = (ArrayList<PhyloTree>)ois.readObject();
        System.out.println("Loading trees trifurcations");
        ArrayList<ArrayList<PhyloTree>> experimentTreesTrifurcations = (ArrayList<ArrayList<PhyloTree>>)ois.readObject();

        System.out.println("################################################");
        System.out.println("NxIndex="+NxIndex);
        System.out.println("pruningIndex="+NxIndex);
        System.out.println("expectedPlacements="+expectedPlacements);
        System.out.println("prunedTrees="+experimentTrees.size());
        for (int i = 0; i < experimentTrees.size(); i++) {
            PhyloTree get = experimentTrees.get(i);
            System.out.println(i+"th tree size test:"+get.getNodeCount());
        }
        System.out.println("################################################");

        //load Dtx
        System.out.println("Loading Dtx matrix...");
        File DtxFile = new File(workDir + File.separator + "Dtx.csv");
        Dtx Dtx=new Dtx(DtxFile);
        //System.out.println(Dtx);

        /// Find all .jplace files in the list of placement results directories
        List<Path> allJplaceFiles = findJplaceFiles(dirs);

        /// Create the set of columns for results.csv based on jplace files
        var paramSet = initializeParamSet(allJplaceFiles);

        // Prepare the results.csv file
        Path csvResult = Paths.get(workDir.getAbsolutePath(),"results.csv");
        BufferedWriter bw= Files.newBufferedWriter(csvResult);

        // Write .csv header
        String parametersJoined = String.join(";", paramSet.keySet());
        bw.append(parametersJoined).append("\n");

        // For each jplace file, calculate node distance
        for (Path currentJPlaceFile : allJplaceFiles) {
            //System.out.println(currentJPlaceFile);
            String jplaceLabel = currentJPlaceFile.getFileName().toString().split("\\.jplace$")[0];
            String[] elts = jplaceLabel.split("_");

            //information related to this placement
            int pruning = Integer.valueOf(elts[0]);
            String software = elts[elts.length - 1];

            System.out.println("--------------------------------------");
            System.out.println("Parsing " + currentJPlaceFile.getFileName().toString());
            System.out.println("software:" + software + " pruning:" + pruning);

            String[] infos = jplaceLabel.split("_");
            System.out.println("Run parameters: " + Arrays.toString(infos));
            TreeMap<Integer, String> paramsValues = new TreeMap<>();
            for (int idx = 1; idx < infos.length - 1; idx++) {
                Matcher m = patternMatch.matcher(infos[idx]);
                if (m.matches()) {
                    String param = m.group(1);
                    if (paramSet.containsKey(param)) {
                        String val = m.group(2);
                        paramsValues.put(paramSet.get(param), val);
                    }
                } else {
                    System.out.println("Error in jplace filename parameters coding, do not matches expected pattern: " + infos[idx]);
                    System.exit(1);
                }
            }
            //paramsValues.put(paramSet.get("file"),jplaceLabel);
            paramsValues.put(paramSet.get("pruning"), Integer.toString(pruning));
            paramsValues.put(paramSet.get("software"), software);

            //System.out.println(paramsValues);

            //load tree and expectedPlacements related to this pruning
            PhyloTree experimentTree = experimentTrees.get(pruning);
            experimentTree.initIndexes();
            ArrayList<Integer> experimentPlacements = expectedPlacements.get(pruning);

            //load jplace content
            JplacerLoader jpl = null;
            if (!software.equals("epang")) {
                jpl = new JplacerLoader(currentJPlaceFile.toFile(), false);
            } else {
                //was false with older versions of epang, but now seems to conserve root correctly,
                //code kept if need to reverse
                jpl = new JplacerLoader(currentJPlaceFile.toFile(), false);
            }


            if (jpl.getTree().getNodeCount() != experimentTree.getNodeCount()) {
                jpl.getTree().displayTree();
                experimentTree.displayTree();
                System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                System.out.println("They do not include the same trees for the same Nx experiment.");
                //return;
                System.exit(1);
            }

            //version of EPA-ng prior to 0.3.4 unroots the input tree
            //this needs to be corrected, posterior versions, no need to correct

            if (jpl.getTree().getNodeCount() != experimentTree.getNodeCount()) {
                System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                System.out.println("They do not include the same trees for the same Nx experiment.");
                System.exit(1);
            }

            //map EPA jplace to experimentTree
            HashMap<Integer, Integer> mapNodes = jpl.getTree().mapNodes(experimentTree);
            //retrieve best placements
            HashMap<String, ArrayList<Placement>> bestPlacements = jpl.getPlacements();

            //iterate on placements
            for (Iterator<String> iterator = bestPlacements.keySet().iterator(); iterator.hasNext(); ) {
                String name = iterator.next();

                int topND = -1;
                double topLwr = -1;
                ArrayList<Integer> nds = new ArrayList<>();
                ArrayList<Double> lwrs = new ArrayList<>();
                double lwrSum = 0.0;
                //iterate on placement branches once to compute (expected_ND)*LWR sum
                for (int pla = 0; pla < bestPlacements.get(name).size(); pla++) {
                    //get the best placement as the nodeId of the phylotree generated
                    //during jplace parsing
                    Integer jplacePhyloTreeNodeId = bestPlacements.get(name).get(pla).getNodeId();
                    double lwr = bestPlacements.get(name).get(pla).getWeightRatio();
                    if (pla == 0) {
                        topLwr = lwr;
                    }
                    //get its equivalent nodeId in the phylotree loaded from the
                    //expected_placements.bin
                    int experimentTreeNodeId = mapNodes.get(jplacePhyloTreeNodeId);
                    //calculate the distance between these 2 nodeIds
                    //i.e. use the DTx matrix
                    int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), experimentTreeNodeId);
                    if (pla == 0) {
                        topND = nodeDistance;
                    }
                    lwrSum += lwr;
                }

                //iterate again to compute eND
                double expectedNodeDistance = 0.0;
                for (int pla = 0; pla < bestPlacements.get(name).size(); pla++) {
                    //get best placement as the nodeId of the phylotree generated
                    //during jplace parsing
                    Integer jplacePhyloTreeNodeId = bestPlacements.get(name).get(pla).getNodeId();
                    double lwr = bestPlacements.get(name).get(pla).getWeightRatio();
                    if (pla == 0) {
                        topLwr = lwr;
                    }
                    //get its equivalent nodeId in the phylotree loaded from the
                    //expected_placements.bin
                    int experimentTreeNodeId = mapNodes.get(jplacePhyloTreeNodeId);
                    //calculate the distance between these 2 nodeIds
                    //i.e. use the DTx and D'Tx matrices
                    int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), experimentTreeNodeId);
                    if (pla == 0) {
                        topND = nodeDistance;
                    }
                    expectedNodeDistance += nodeDistance * lwr / lwrSum;
                }

                //get coordinates of placed read in original alignment (before pruning)
                String[] readInfos = name.split("_");
                //System.out.println("Read name: " + name);
                long readStart = 0;
                long readEnd = 0;
                try {
                    readStart = Long.decode(readInfos[readInfos.length - 2]);
                } catch (NumberFormatException ex) {
                    ex.printStackTrace();
                }
                try {
                    readEnd = Long.decode(readInfos[readInfos.length - 1]);
                } catch (NumberFormatException ex) {
                    ex.printStackTrace();
                }

                //add information only in specific columns
                paramsValues.put(paramSet.get("query"), name);
                paramsValues.put(paramSet.get("rstart"), Long.toString(readStart));
                paramsValues.put(paramSet.get("rend"), Long.toString(readEnd));
                paramsValues.put(paramSet.get("nd"), Integer.toString(topND));
                paramsValues.put(paramSet.get("e_nd"), Double.toString(expectedNodeDistance));

                System.out.println(topND);


                //now build output string
                for (int column = 0; column < paramSet.keySet().size(); column++) {
                    if (column > 0) {
                        bw.append(";");
                    }
                    if (paramsValues.containsKey(column)) {
                        bw.append(paramsValues.get(column));
                    }
                }
                bw.append("\n");
            }


        }

        bw.close();
    }

    public static void main(String[] args) {
        
        //System.setProperty("debug.verbose","1");

        try {
            
            System.out.println("ARGS: workDir [list_of_tested_software_directories,comma-separated]");
            System.out.println("example: /path/to/pewo_workdir EPANG,RAPPAS,PPLACER");
            System.out.println(System.getProperty("java.class.path"));

            if(args.length>0) {
                File workingDirectory = new File(args[0]);
                ArrayList<File> placementSoftwareDirs = new ArrayList<>(20);

                System.out.println("workDir: " + workingDirectory);
                if (! workingDirectory.exists()) {
                    System.out.println("workdir does noot exists");
                    System.exit(1);
                }

                String[] to_test = args[1].split(",");
                for (String s : to_test) {
                    File expected = new File(workingDirectory + File.separator + s);

                    if (!expected.exists()) {
                        System.out.println("Directory does not exists: " + expected.getAbsolutePath());
                        System.exit(1);
                    } else {
                        placementSoftwareDirs.add(expected);
                    }
                }

                // launch
                DistanceGenerator dg = new DistanceGenerator(workingDirectory, placementSoftwareDirs);
                dg.computeDistances();

            } else {
                System.out.println("No arguments passed!");
                System.exit(1);
            }

        } catch (IOException | ClassNotFoundException ex) {
            Logger.getLogger(DistanceGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}
