import java.util.*;

public class MovingAverage {

    public static SortedMap<Double,Double> averagedMap(SortedMap<Double,Double> hashMap, int windowSize){
        double sum;
        double[] values = hashMap.values().stream().mapToDouble(Double::doubleValue).toArray();

        SortedMap<Double, Double> output = new TreeMap<>();
        int position = -1;
        Iterator<Double> forward = hashMap.keySet().iterator();
        while(forward.hasNext()){
            int count = 1;
            ++position;
            sum = values[position];
            for (int i = 1; i <= windowSize; i++) {
                if(position-i >= 0){
                    sum += values[position-i];
                    ++count;
                }
                if(position+i < values.length){
                    sum += values[position+i];
                    ++count;
                }
            }
            output.put(forward.next(), (sum/count));
        }
        return output;
    }

    //Would love to average based on a method like that below, but the computation time is several orders of magnitude greater
/*
        public Map<Double,Double> averagedMap(double range) {
        Map<Double, Double> output = Collections.synchronizedMap(new HashMap<>());
        hashMap.forEach((key,value) -> {
            Double fromKey = key - (range);
            Double toKey = key + (range);
            double newkey = hashMap.subMap(fromKey, toKey).keySet().stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double newvalue = hashMap.subMap(fromKey, toKey).values().stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            synchronized (output){output.put(newkey, newvalue);}
        });
        return output;
    }*/
}
