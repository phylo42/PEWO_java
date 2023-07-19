
import dtx.Dtx;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
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
    final Pattern classicParameterPattern = Pattern.compile("([a-z]+)([0-9A-z\\-\\.]+)");
    final Pattern localGroupPattern = Pattern.compile("^([a-z0-9]+[A-Z0-9\\.]+)(-[a-z0-9]+[A-Z0-9\\.]+)*$");

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

    /// An auxilary class used to parse parameter key-value pairs in .jplace file names
    static class Parameter {
        public String name;
        public String value;

        public Parameter(String n, String v) {
            name = n;
            value = v;
        }
    };

    /// Parses a group of local parameters. Returns a list of key-value pairs with
    /// local parameter names resolved using the value of the first parameter.
    /// Example: rgenART-coILLUMINA-plMSV3 is resolved to
    /// [ (art-co, ILLUMINA), (art-pl, MSV3) ]
    List<Parameter> parseLocalGroup(String parameterGroup) {
        ArrayList<Parameter> parameters = new ArrayList<>();

        String[] keyValuePairs = parameterGroup.split("-");

        /// Form the group prefix, i.e. take the value of the first
        /// local parameter as the prefix for the whole group of local paramters.
        /// E.g. the group "rgenART-coILLUMINA-plHS25" will have
        /// parameters art-co, art-pl (i.e. group prefix "art")
        String groupPrefix = keyValuePairs[0].split("(?<=\\p{Lower})(?=\\p{Upper})")[1].toLowerCase();

        /// Iterate over local parameters but the first one
        for (int j = 1; j < keyValuePairs.length; j++) {
            String keyValuePair = keyValuePairs[j];
            String[] keyAndValue = keyValuePair.split("(?<=\\p{Lower})(?=\\p{Upper})");
            String key = keyAndValue[0];
            String value = keyAndValue[1];

            /// Resolve the local parameter name to the global space:
            /// e.g. art-co, art-pl
            String resolvedName = groupPrefix + "-" + key;
            parameters.add(new Parameter(resolvedName, value));
        }
        return parameters;
    }

    /// Updates the paramSet with the columns specific to tested
    /// software based on jplace file names
    void initializeUncommonColumns(List<Path> jplaceFiles, int columnCounter) {


        for (Path currentJplaceFile : jplaceFiles) {
            String filename = currentJplaceFile.toFile().getName();

            /// Split the filename by _ into parameter groups.
            /// The classic PEWO convention is that all parameters
            /// are just separated by the underscore: e.g. 0_r150_ms6_sb3_mp40_pplacer.jplace
            /// The new convention assumes that underscores split the name
            /// into "scopes". Every scope can be either a global key-value parameter (classic)
            /// or a dash-concatenated sequences of local key-value parameters (new).
            String[] parameterGroups = filename.split("_");

            /// 1st element is pruning id, last element is "software.jplace"
            for (int i = 1; i < parameterGroups.length - 1; i++) {
                String parameterGroup = parameterGroups[i];

                /// New convention
                if (parameterGroup.contains("-")) {
                    if (!localGroupPattern.matcher(parameterGroup).matches()) {
                        throw new RuntimeException("Invalid format: " + parameterGroup);
                    }

                    List<Parameter> localParams = parseLocalGroup(parameterGroup);
                    for (Parameter parameter : localParams) {
                        if (!paramSet.containsKey(parameter.name)) {
                            paramSet.put(parameter.name, ++columnCounter);
                        }
                    }


                }
                /// Classic convention
                else
                {
                    Matcher m = classicParameterPattern.matcher(parameterGroup);
                    if (m.matches()) {
                        String param = m.group(1);
                        if (!paramSet.containsKey(param)) {
                            paramSet.put(m.group(1), ++columnCounter);
                        }
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

    static class ReadParameters
    {
        public long readStart = 0;
        public long readEnd = 0;
    };

    ReadParameters parseReadParameters(String queryName) {
        ReadParameters parameters = new ReadParameters();

        /// Get coordinates of placed read in original alignment (before pruning)
        String[] readInfo = queryName.split("_");
        try {
            parameters.readStart = Long.decode(readInfo[readInfo.length - 2]);
            parameters.readEnd = Long.decode(readInfo[readInfo.length - 1]);
        } catch (NumberFormatException | ArrayIndexOutOfBoundsException ex) {
            /// Reads generated by programs like ART do not have neither start
            /// nor end. An exception is expected
        }
        return parameters;
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
        LinkedHashMap<String,Integer> paramSet = initializeParamSet(allJplaceFiles);

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
                String parameterGroup = infos[idx];

                Matcher m = classicParameterPattern.matcher(parameterGroup);
                if (m.matches()) {
                    /// New convention with local parameters
                    if (parameterGroup.contains("-")) {
                        List<Parameter> localParameters = parseLocalGroup(parameterGroup);

                        for (Parameter parameter : localParameters) {
                            if (paramSet.containsKey(parameter.name)) {
                                paramsValues.put(paramSet.get(parameter.name), parameter.value);
                            }
                        }
                    }
                    /// Classic convention
                    else {
                        String param = m.group(1);
                        String val = m.group(2);

                        if (paramSet.containsKey(param)) {
                            paramsValues.put(paramSet.get(param), val);
                        }
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

            // load tree and expectedPlacements related to this pruning
            PhyloTree experimentTree = experimentTrees.get(pruning);
            experimentTree.initIndexes();
            ArrayList<Integer> experimentPlacements = expectedPlacements.get(pruning);

            // load jplace content
            JplacerLoader jpl = null;
            if (!software.equals("epang")) {
                jpl = new JplacerLoader(currentJPlaceFile.toFile(), false);
            } else {
                // was false with older versions of epang, but now seems to conserve root correctly,
                // code kept if needed to reverse
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

            // version of EPA-ng prior to 0.3.4 unroots the input tree
            // this needs to be corrected, posterior versions, no need to correct

            if (jpl.getTree().getNodeCount() != experimentTree.getNodeCount()) {
                System.out.println("Something is wrong between the JPlace and expected_placements.bin trees.");
                System.out.println("They do not include the same trees for the same Nx experiment.");
                System.exit(1);
            }

            // map EPA jplace to experimentTree
            HashMap<Integer, Integer> mapNodes = jpl.getTree().mapNodes(experimentTree);
            // retrieve best placements
            HashMap<String, ArrayList<Placement>> bestPlacements = jpl.getPlacements();

            //iterate on placements
            for (String queryName : bestPlacements.keySet()) {
                paramsValues.put(paramSet.get("query"), queryName);

                int topND = -1;
                double topLwr = -1;
                double lwrSum = 0.0;
                // iterate on placement branches once to compute (expected_ND)*LWR sum
                for (int pla = 0; pla < bestPlacements.get(queryName).size(); pla++) {
                    // get the best placement as the nodeId of the phylotree generated
                    // during jplace parsing
                    Integer jplacePhyloTreeNodeId = bestPlacements.get(queryName).get(pla).getNodeId();
                    double lwr = bestPlacements.get(queryName).get(pla).getWeightRatio();
                    if (pla == 0) {
                        topLwr = lwr;
                    }
                    // get its equivalent nodeId in the phylotree loaded from the
                    // expected_placements.bin
                    int experimentTreeNodeId = mapNodes.get(jplacePhyloTreeNodeId);
                    // calculate the distance between these 2 nodeIds
                    // i.e. use the DTx matrix
                    int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), experimentTreeNodeId);
                    if (pla == 0) {
                        topND = nodeDistance;
                    }
                    lwrSum += lwr;
                }

                // iterate again to compute eND
                double expectedNodeDistance = 0.0;
                for (int pla = 0; pla < bestPlacements.get(queryName).size(); pla++) {
                    // get the best placement as the nodeId of the phylotree generated
                    // during jplace parsing
                    Integer jplacePhyloTreeNodeId = bestPlacements.get(queryName).get(pla).getNodeId();
                    double lwr = bestPlacements.get(queryName).get(pla).getWeightRatio();
                    if (pla == 0) {
                        topLwr = lwr;
                    }
                    // get its equivalent nodeId in the phylotree loaded from the
                    // expected_placements.bin
                    int experimentTreeNodeId = mapNodes.get(jplacePhyloTreeNodeId);

                    // calculate the distance between these 2 nodeIds
                    // i.e. use the DTx and D'Tx matrices
                    int nodeDistance = Dtx.getNodeDistance(pruningIndex.get(pruning), experimentTreeNodeId);
                    if (pla == 0) {
                        topND = nodeDistance;
                    }
                    expectedNodeDistance += nodeDistance * lwr / lwrSum;
                }
                paramsValues.put(paramSet.get("nd"), Integer.toString(topND));
                paramsValues.put(paramSet.get("e_nd"), Double.toString(expectedNodeDistance));

                // Parse read_start, read_end if they present in the query name
                ReadParameters readParameters = parseReadParameters(queryName);
                paramsValues.put(paramSet.get("rstart"), Long.toString(readParameters.readStart));
                paramsValues.put(paramSet.get("rend"), Long.toString(readParameters.readEnd));

                // Make the output string
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
