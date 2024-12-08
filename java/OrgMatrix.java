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

import edu.uw.bothell.css.dsl.MASS.MASSBase;
import edu.uw.bothell.css.dsl.MASS.MASS;
import edu.uw.bothell.css.dsl.MASS.VertexPlace;
import edu.uw.bothell.css.dsl.MASS.Place;

import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.Collectors;

public class OrgMatrix extends Place {
    private String keyName = "";
    private HashMap<String, Double> map = new HashMap<>();
    private double keyVal = 0;
    private boolean status = true;

    // function identifiers
    public static final int init_ = 0;
    public static final int findMin_ = 1;
    public static final int updateMap_ = 2;
    public static final int updateNewKey_ = 3;

    public Object callMethod(int functionId, Object argument) {
        switch (functionId) {
            case init_:
                return init(argument);
            case updateMap_:
                return updateMap(argument);
            case findMin_:
                return findMin();
            case updateNewKey_:
                return updateNewKey(argument);
        }
        return null;
    }

    private double getPrevVal(String key) {
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

    private String formatToSixDecimals(double number) {
        return String.format("%.6f", number);
    }

    private Object[] updateMap(Object arugment) {
        if(!status) {
            return new String[] {"invalid"};
        }
        String[] pair = ((String)arugment).split(" ");

        String key1 = pair[0];
        String key2 = pair[1];
        double prevKey1Val = getPrevVal(key1);
        double prevKey2Val = getPrevVal(key2);
        double val = Double.parseDouble(pair[2]) / 2;

        if(key2.equals(keyName) || key1.equals(keyName)) {
            status = false;

            return new String[] {"invalid"};
        }

        String newKey = "(" + key1 + ":" + formatToSixDecimals(val - prevKey1Val) + "," + key2 + ":" + formatToSixDecimals(val - prevKey2Val) + ")";
        
        int sizeOfKey1 = key1.split(",").length;
        int sizeOfKey2 = key2.split(",").length;
        double localValOfKey1 = map.get(key1);
        double localValOfKey2 = map.get(key2);

        double newVal = (localValOfKey1 * sizeOfKey1 + localValOfKey2 * sizeOfKey2) / (sizeOfKey1 + sizeOfKey2);
        map.remove(key1);
        map.remove(key2);
        map.put(newKey, newVal);

        return new String[] { keyName + " " + String.valueOf(newVal) };
    }


    // 0: prevKey1Name, newKeyName, updated value
    private Object updateNewKey(Object argument) {
        String[] args = (String[])argument;
        
        String key = args[0];

        if(!key.equals(keyName)) {
            return "";
        }

        status = true;
        String newKey = args[1];
        keyName = newKey;
        // System.out.println("newKey: " + newKey);
        map.clear();

        for(int i = 2; i < args.length; i++) {
            String[] pair = args[i].split(" ");
            map.put(pair[0], Double.parseDouble(pair[1]));
        }

        return "";
    }

    private Object findMin() {
        if(!status) {
            return new String[] {"", "", "-1"};
        }

        double min = Double.MAX_VALUE;
        String localKeyName = "";
        for(String localKey : map.keySet()) {
            double localVal = map.get(localKey);

            if(localVal < min) {
                min = localVal;
                localKeyName = localKey;
            }
        }

        if(keyName.compareTo(localKeyName) == 1) {
            return new String[] { localKeyName, keyName, String.valueOf(min) };
        }
        
        return new String[] { keyName, localKeyName, String.valueOf(min) };
    }

    public Object init(Object argument) {

        int rowIndex = this.getIndex()[0];
        String text = (String) argument;

        String[] token = text.trim().split("\\s+");
        keyName = token[0].trim();

        for(int i = 1; i < token.length; i++) {
            String localKey = "Org" + i;
            if(localKey.equals(keyName)){
                continue;
            }

            map.put(localKey, Double.parseDouble(token[i]));
        }
        return "";
    }

    /**
     * Is the default constructor.
     */
    public OrgMatrix() {
        // super();
    }

    public OrgMatrix(Object arg) {
        // super(arg);
    }
}
