import java.io.IOException;
import org.apache.hadoop.conf.Configuration;
import java.io.InputStreamReader;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.fs.FileSystem;
import java.io.BufferedReader;
import java.io.FileReader;
import java.net.URI;
import org.apache.hadoop.filecache.DistributedCache;
import java.util.*;

public class UPGMAMapReduce {

    public static class DistanceMapper extends Mapper<Object, Text, Text, Text> {
        @Override
        public void map(Object key, Text value, Context context) throws IOException, InterruptedException {
            String[] tokens = value.toString().trim().split("\\s+");
            if (tokens.length <= 2) {
                return;
            }

            boolean hasMin = false;
            String minKey = "";
            double minVal = Double.MAX_VALUE;
            String mainKey = tokens[0].trim();

            for (int i = 1; i < tokens.length; i++) {
                String[] pair = tokens[i].trim().split(":");
                String localKey = pair[0].trim();

                if (mainKey.equals(localKey)) {
                    continue;
                }

                double localVal = Double.parseDouble(pair[1].trim());
                if (localVal <= minVal) {
                    hasMin = true;
                    minKey = localKey;
                    minVal = localVal;
                }
            }

            if (hasMin) {
                String pair = (mainKey.compareTo(minKey) == 1) ?
                    (mainKey + " " + minKey + " " + minVal) :
                    (minKey + " " + mainKey + " " + minVal);
                context.write(new Text("min"), new Text(pair));
            }
        }
    }

    public static class MinimumReducer extends Reducer<Text, Text, Text, Text> {
        @Override
        public void reduce(Text key, Iterable<Text> values, Context context) throws IOException, InterruptedException {
            double minVal = Double.MAX_VALUE;
            String minPair = "";
            boolean hasMin = false;

            // 'globalMin' is not defined anywhere in this code snippet.
            // Assuming it should be defined or calculated before use.
            // For now, let's assume globalMin is meant to be minVal or is a known variable.
            double globalMin = 0; // Placeholder: Adjust according to actual logic

            for (Text val : values) {
                String[] tokens = val.toString().split("\\s+");
                double localVal = Double.parseDouble(tokens[2].trim());
                if (localVal <= minVal) {
                    hasMin = true;
                    minVal = localVal;
                    minPair = tokens[0].trim() + " " + tokens[1].trim() + " " + minVal;
                }
            }

            if (hasMin) {
                context.write(new Text(minPair), new Text(String.valueOf(globalMin / 2)));
            }
        }
    }

    public static class UpdateMapper extends Mapper<Object, Text, Text, Text> {
        private double keyValue = 0;
        private String key1 = "";
        private String key2 = "";
        private double prevKey1Val = 0;
        private double prevKey2Val = 0;
        private String newKey = "";
        private int nOfKey1 = 1;
        private int nOfKey2 = 1;

        @Override
        protected void setup(Context context) throws IOException, InterruptedException {
            Path filePath = new Path("hdfs:///user/ypang5/output/part-r-00000");
            FileSystem fs = FileSystem.get(context.getConfiguration());

            if (fs == null) {
                throw new IOException("FileSystem is null. Check Hadoop configuration.");
            }
            if (!fs.exists(filePath)) {
                throw new IOException("File does not exist: " + filePath);
            }

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(fs.open(filePath)))) {
                String line = reader.readLine();
                if (line != null && !line.isEmpty()) {
                    String[] tokens = line.split("\\s+");
                    key1 = tokens[0].trim();
                    key2 = tokens[1].trim();
                    prevKey1Val = getPrevKeyVal(key1);
                    prevKey2Val = getPrevKeyVal(key2);
                    keyValue = Double.parseDouble(tokens[2]);
                    nOfKey1 = key1.split(",").length;
                    nOfKey2 = key2.split(",").length;
                    newKey = "(" + key1 + "|" + formatToSixDecimals((keyValue / 2 - prevKey1Val)) + "," +
                             key2 + "|" + formatToSixDecimals((keyValue / 2 - prevKey2Val)) + ")";
                }
            }
        }

        private double getPrevKeyVal(String key) {
            if (key.contains("|")) {
                String[] tokens = key.substring(1, key.length() - 1).split(",");
                String[] pairs1 = tokens[0].split("\\|");
                String[] pairs2 = tokens[1].split("\\|");
                return Math.max(Double.parseDouble(pairs1[pairs1.length - 1]), Double.parseDouble(pairs2[pairs2.length - 1]));
            }
            return 0;
        }

        private String formatToSixDecimals(double number) {
            return String.format("%.6f", number);
        }

        @Override
        public void map(Object key, Text value, Context context) throws IOException, InterruptedException {
            String[] tokens = value.toString().trim().split("\\s+");
            String mainKey = tokens[0].trim();

            if (mainKey.equals(key1) || mainKey.equals(key2)) {
                context.write(new Text("newKey"), new Text(""));
                return; // No output for these keys
            }

            HashMap<String, Double> map = new HashMap<>();

            for (int i = 1; i < tokens.length; i++) {
                String[] pair = tokens[i].trim().split(":");
                map.put(pair[0].trim(), Double.parseDouble(pair[1].trim()));
            }

            double valKey1 = map.getOrDefault(key1, 1000000.0);
            double valKey2 = map.getOrDefault(key2, 1000000.0);

            double newDistance = 1000000.0;
            if (valKey1 != 1000000 && valKey2 != 1000000) {
                newDistance = (valKey1 * nOfKey1 + valKey2 * nOfKey2) / (nOfKey1 + nOfKey2);
            }

            map.put(newKey, newDistance);
            context.write(new Text("newKey"), new Text(mainKey + ":" + formatToSixDecimals(newDistance)));

            StringBuilder sb = new StringBuilder();
            sb.append(mainKey);
            sb.append(" ");
            for (String localKey : map.keySet()) {
                if (localKey.equals(key1) || localKey.equals(key2)) {
                    continue;
                }
                double localVal = map.get(localKey);
                sb.append(localKey).append(":").append(localVal).append(" ");
            }

            context.write(new Text("existingKey"), new Text(sb.toString()));
        }
    }

    public static class UpdateReducer extends Reducer<Text, Text, Text, Text> {
        private double keyValue = 0;
        private String key1 = "";
        private String key2 = "";
        private double prevKey1Val = 0;
        private double prevKey2Val = 0;
        private String newKey = "";
        private int nOfKey1 = 1;
        private int nOfKey2 = 1;

        @Override
        protected void setup(Context context) throws IOException, InterruptedException {
            Path filePath = new Path("hdfs:///user/ypang5/output/part-r-00000");
            FileSystem fs = FileSystem.get(context.getConfiguration());

            if (fs == null) {
                throw new IOException("FileSystem is null. Check Hadoop configuration.");
            }
            if (!fs.exists(filePath)) {
                throw new IOException("File does not exist: " + filePath);
            }

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(fs.open(filePath)))) {
                String line = reader.readLine();
                if (line != null && !line.isEmpty()) {
                    String[] tokens = line.split("\\s+");
                    key1 = tokens[0].trim();
                    key2 = tokens[1].trim();
                    prevKey1Val = getPrevKeyVal(key1);
                    prevKey2Val = getPrevKeyVal(key2);
                    keyValue = Double.parseDouble(tokens[2]);
                    nOfKey1 = key1.split(",").length;
                    nOfKey2 = key2.split(",").length;
                    newKey = "(" + key1 + "|" + formatToSixDecimals((keyValue / 2 - prevKey1Val)) + "," +
                             key2 + "|" + formatToSixDecimals((keyValue / 2 - prevKey2Val)) + ")";
                }
            }
        }

        private double getPrevKeyVal(String key) {
            if (key.contains("|")) {
                String[] tokens = key.substring(1, key.length() - 1).split(",");
                String[] pairs1 = tokens[0].split("\\|");
                String[] pairs2 = tokens[1].split("\\|");
                return Math.max(Double.parseDouble(pairs1[pairs1.length - 1]), Double.parseDouble(pairs2[pairs2.length - 1]));
            }
            return 0;
        }

        private String formatToSixDecimals(double number) {
            return String.format("%.6f", number);
        }

        @Override
        public void reduce(Text key, Iterable<Text> values, Context context) throws IOException, InterruptedException {
            if (key.toString().equals("newKey")) {
                StringBuilder sb = new StringBuilder();
                sb.append(newKey).append(" ");
                for (Text val : values) {
                    sb.append(val.toString()).append(" ");
                }
                sb.append(newKey).append(":0");
                context.write(null, new Text(sb.toString()));
                return;
            }

            for (Text val : values) {
                context.write(null, new Text(val));
            }
        }
    }

public static void main(String[] args) throws Exception {
    Configuration conf = new Configuration();
    FileSystem fs = FileSystem.get(conf);
    Path outputPath = new Path(args[1]);
    Path inputPath = new Path(args[0]);
    Path tempPath = new Path(args[2]);

    // Initial cleanup if needed
    if (fs.exists(outputPath)) {
        fs.delete(outputPath, true);
    }
    if (fs.exists(tempPath)) {
        fs.delete(tempPath, true);
    }
   
    int step = 0;
    long jobStartTime = System.currentTimeMillis(); // Record start time

    while (true) {
        // First job
        Job job = new Job(conf, "WordCount");
        job.setJarByClass(WordCount.class);
        job.setMapperClass(DistanceMapper.class);
        job.setReducerClass(MinimumReducer.class);
        job.setOutputKeyClass(Text.class);
        job.setOutputValueClass(Text.class);

        FileInputFormat.addInputPath(job, inputPath);
        FileOutputFormat.setOutputPath(job, outputPath);

        if (!job.waitForCompletion(true)) {
            System.err.println("Job 1 failed.");
            System.exit(1);
        }

        // Second job
        Job job2 = new Job(conf, "updateTable");
        job2.setJarByClass(WordCount.class);
        job2.setMapperClass(UpdateMapper.class);
        job2.setReducerClass(UpdateReducer.class);
        job2.setOutputKeyClass(Text.class);
        job2.setOutputValueClass(Text.class);

        FileInputFormat.addInputPath(job2, inputPath);
        FileOutputFormat.setOutputPath(job2, tempPath);

        if (!job2.waitForCompletion(true)) {
            System.err.println("Job 2 failed.");
            System.exit(1);
        }
       
       
        System.out.println("current step: " + step);
        step++;
       
        if (fs.exists(inputPath)) {
            fs.delete(inputPath, true);
        }

            // Delete old input and output so next iteration is clean
        if (fs.exists(outputPath)) {
            fs.delete(outputPath, true);
        }

        fs.rename(tempPath, inputPath);

        // Check line count in tempPath
        Path resultFile = new Path(inputPath, "part-r-00000");
        int lineCount = 0;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(fs.open(resultFile)))) {
            while (br.readLine() != null) {
                lineCount++;
                if (lineCount > 1) {
                    break; // We only need to know if it's more than one line
                }
            }
        }

        if (lineCount <= 1) {
            // Only one line - finish the program
            System.out.println("Only one line in the tempPath file. Finishing job...");
           
           
      long jobEndTime = System.currentTimeMillis();
        long elapsedTimeMillis = jobEndTime - jobStartTime;
        double elapsedTimeSeconds = elapsedTimeMillis;
        System.out.println("This iteration took: " + elapsedTimeSeconds + " seconds.");
            System.exit(0);
        }
    }
   

}

}