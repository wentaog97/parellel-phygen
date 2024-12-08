/*

 	MASS Java Software License
	© 2012-2021 University of Washington

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in
	all copies or substantial portions of the Software.

	The following acknowledgment shall be used where appropriate in publications, presentations, etc.:      

	© 2012-2021 University of Washington. MASS was developed by Computing and Software Systems at University of 
	Washington Bothell.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
	THE SOFTWARE.

*/

package edu.uw.bothell.css.dsl.appl.graphs.TriangleCounting;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

import edu.uw.bothell.css.dsl.MASS.Agents;
import edu.uw.bothell.css.dsl.MASS.GraphPlaces;
import edu.uw.bothell.css.dsl.MASS.MASS;
import edu.uw.bothell.css.dsl.MASS.MASSBase;
import edu.uw.bothell.css.dsl.MASS.graph.CytoscapeListener;
import edu.uw.bothell.css.dsl.MASS.graph.Graph;
import edu.uw.bothell.css.dsl.MASS.graph.MASSListener;
import edu.uw.bothell.css.dsl.MASS.graph.transport.VertexModel;
import edu.uw.bothell.css.dsl.MASS.logging.LogLevel;
import edu.uw.bothell.css.dsl.appl.graphs.TriangleCounting.OrgMatrix;
import edu.uw.bothell.css.dsl.MASS.Places;


/**
 * CountTriangles.java - Counts and enumerates all triangles in a given graph
 *
 * @author Munehiro Fukuda <mfukuda@uw.edu>
 * @version 1.0
 * @since October 23, 2017
 */
public class UPGMAMASS {
    private static final boolean CYTOSCAPE_LISTENER = false;
    private static int sequentialId = 2;

    /**
     * Counts and enumerates all triangles in a given graph
     *
     * @param args TBD
     */
    public static void main(String[] args) {
        // TODO: This program _appears_ to have 1 invariant:
        //     1. The number of vertices must be correct

        // Read and validate input parameters
        if (args.length != 1) {
            System.err.println("Usage: java -jar TriangleCounting.jar <dsl_graph_file>");
            System.exit(-1);
        }
        
        String filePath = args[0];

        System.out.println("filePath: " + filePath);

        MASS.setLoggingLevel(LogLevel.DEBUG);
        MASS.init(100000);

        System.out.println("Begin data space generation");
        
        ArrayList<String> inputList = loadDataFromFile(filePath);
        String[] globalData = inputList.toArray(new String[0]);

        for(String item: globalData) {
            System.out.println(item);
        }

        long begin = System.currentTimeMillis();
        // int MAX = Integer.MAX_VALUE;

        int numRows = globalData.length;

        Places places = new Places(1, OrgMatrix.class.getName(), (Object)new Integer(0), numRows);

        Object[] res = places.callAll(OrgMatrix.init_, globalData);
        
        // boolean execu = false;

        int iteration = 1;
        while(true) {
        // GraphPlaces network = new GraphPlaces(1, NodeGraphMASS.class.getName());
            try {
                Object[] minItems = places.callAll(OrgMatrix.findMin_, globalData);
                // System.out.println("findMin_ Length: " + minItems.length);

                double min = Double.MAX_VALUE;
                String[] minPair = new String[] {"", "", ""};     
                for(int i = 0; i < minItems.length; i++) {
                    String[] item = (String[])minItems[i];
                    double val = Double.parseDouble(item[2]);
                    if(val < 0) {
                        continue;
                    }

                    if(min > val) {
                        min = val;
                        minPair = item;
                    }
                }

                if(minPair[0].equals("")){
                    System.out.println("incorrect result");
                    break;
                }


                if(minPair[0].split(",").length + minPair[1].split(",").length == numRows) {
                    System.out.println("new wick tree: " + generateNewKey(minPair[0], minPair[1], minPair[2]));
                    break;
                }

                String minArgs = minPair[0] + " " + minPair[1] + " " + minPair[2];

                System.out.println("min val:" + minPair[0] + "; " + minPair[1] + "; " + minPair[2]);

                // take a look debug here
                Object[] updatedItems = places.callAll(OrgMatrix.updateMap_, new String[] {minArgs});
                places.callAll(OrgMatrix.updateMap_, new String[] {minArgs});

                ArrayList<String> newKeyItems = new ArrayList<>();
                newKeyItems.add(minPair[0]);
                newKeyItems.add(generateNewKey(minPair[0], minPair[1], minPair[2]));
                String[] items = (String[])updatedItems;
                for(int i = 0; i < items.length; i++) {
                    String item = items[i];

                    if(item.equals("invalid")) {
                        continue;
                    }
                    
                    newKeyItems.add(item);
                }
                places.callAll(OrgMatrix.updateNewKey_, newKeyItems);
                System.out.println("iteration: " + iteration++);
            } catch(Exception e) {
                System.err.println("Error reading in graph file: " + e);
                System.exit(3);
            }
        }

        long end = System.currentTimeMillis();
        
        System.out.println("Import complete\nImport time: " + (end - begin) / 1000.0 + " s");

        System.out.println("Finish MASS");

        MASS.finish();
    }

    private static String generateNewKey(String key1, String key2, String valStr) {
        double prevKey1Val = getPrevVal(key1);
        double prevKey2Val = getPrevVal(key2);
        double val = Double.parseDouble(valStr) / 2;

        return "(" + key1 + ":" + formatToSixDecimals(val - prevKey1Val) + "," + key2 + ":" + formatToSixDecimals(val - prevKey2Val) + ")";
    }

    private static double getPrevVal(String key) {
        if(!key.contains("(")) {
            return 0;
        }

        key = key.substring(1, key.length() - 1);
        String[] tokens = key.split(",");
        String preKey1Str = tokens[0].split(":")[1];
        double preKey1 = Double.parseDouble(preKey1Str);

        String preKey2Val = tokens[1].split(":")[1];
        double preKey2 = Double.parseDouble(preKey2Val);

        return Math.max(preKey1, preKey2);
    }

    private static String formatToSixDecimals(double number) {
        return String.format("%.6f", number);
    }

    public static ArrayList<String> loadDataFromFile(String filename) {
        ArrayList<String> rows = new ArrayList<>();

        try(BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;

            while((line = br.readLine()) != null) {
                // String tokens = line.trim();
                String row = line.trim();
                String updated = row.replace("N/A", "1000000");

                if(!row.equals("")) {
                    rows.add(updated);
                }
            }
        } catch (IOException e) {
            System.err.println("Error reading file: " + e.getMessage());
            e.printStackTrace();
        }

        return rows;
    }
}