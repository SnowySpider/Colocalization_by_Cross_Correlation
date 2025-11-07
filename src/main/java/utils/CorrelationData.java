/*-
 * #%L
 * Scijava plugin for spatial correlation
 * %%
 * Copyright (C) 2019 - 2025 Andrew McCall, University at Buffalo
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package utils;

import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.analysis.function.Gaussian;

import java.util.Iterator;
import java.util.SortedMap;

public class CorrelationData {

    public SortedMap<Double, Double> oCorrMap;
    public SortedMap<Double, Double> sCorrMap;

    public nGaussian gaussians = null;

    //Triplets ordered in: Normalization (Height), Mean, Sigma
    protected double[] gaussFitParameters;

    protected Double[] confidence;
    protected Double rSquared;

    public final int curveCount;

    public double getGaussianNorm(int index){return gaussFitParameters[(index*3)];}
    public double getGaussianMean(int index){return gaussFitParameters[(index*3)+1];}
    public double getGaussianSD(int index){return gaussFitParameters[(index*3)+2];}
    public double getGaussianPeakHeight(int index){return gaussians.getGaussian(index).value(getGaussianMean(index));}
    public boolean hasConfidence(){return confidence != null;}
    public double getConfidence(int index){return confidence[index];}
    public double getRSquared(){return rSquared;}
    public boolean hasGaussians(){return gaussians.getCount() > 0;}

    public CorrelationData(int curveFitCount){
        curveCount = curveFitCount;
    }

    public CorrelationData(double... gaussFitParameters){
        this.gaussFitParameters = gaussFitParameters.clone();
        curveCount = gaussFitParameters.length/3;
        gaussians = new nGaussian(gaussFitParameters);
    }

    public void fitGaussianCurve() {
        gaussFitParameters = CurveFit(sCorrMap);

        gaussians = new nGaussian(gaussFitParameters);

        rSquared = calcRsquared();

        if (oCorrMap != null){
            confidence = new Double[curveCount];
            for (int i = 0; i < curveCount; i++) {
                confidence[i] = (areaUnderCurve(new Gaussian(getGaussianNorm(i), getGaussianMean(i),getGaussianSD(i)), sCorrMap, getGaussianMean(i), getGaussianSD(i)) / areaUnderCurve(oCorrMap, getGaussianMean(i), getGaussianSD(i)));
            }
        }
        for (int i = 0; i < curveCount; i++) {
            if(gaussFitParameters[i+1] == -1)
                throw new NullPointerException("Could not fit Gaussian curve to data");
        }
    }

    private double[] CurveFit(SortedMap<Double, Double> inputMap) {
        double maxLoc = 0;
        double max = 0;
        WeightedObservedPoints obs = new WeightedObservedPoints();

        /*First need to determine the maximum value in order to set the weights for the fitting, and determine its
         * location for instances where the mean is close to zero (in order to mirror the data, this has to be done
         * for a good fit)
         */
        for (Double d : inputMap.keySet()) {
            if (inputMap.get(d) > max) {
                maxLoc = d;
                max = inputMap.get(d);
            }
        }

        /* this next loop adds values below zero that mirror values equidistant from the opposite side of the peak value (max at maxLoc).
         * This is done for fits where the means are near zero, as this data is zero-bounded. Not mirroring the data results
         * in very poor fits for such values. We can't simply mirror across 0 as this will create a double-peak
         * for any data where the peak is near but not at zero.
         * It would be preferable to fit the data using a truncated gaussian fitter, but I could not find any available
         * java class that performs such a fit and my own attempts were unsuccessful.
         */
        if(maxLoc == inputMap.firstKey()){
            maxLoc = 0.0;
        }
        double finalMaxLoc = maxLoc;

        inputMap.forEach((key,value) -> {
            obs.add(key, value);
            if (key > 2 * finalMaxLoc) {
                obs.add(((2 * finalMaxLoc) - key), value);
            }
        });

        double [] output = null;
        nGaussianCurveFitter curveFitter = nGaussianCurveFitter.create();
        try{
            output = curveFitter.withCount(curveCount).withMaxIterations(100).fit(obs.toList());
        }
        catch(Exception ignored){}

        /*
         * Have to check if the curve was fit to a single noise spike, something that came up quite a bit during initial testing.
         * If not fit to a noise spike, values are returned with no further processing, if it is, the data is averaged
         * with nearest neighbors and another fit is attempted. This usually only needs a single round of averaging.
         */

        Iterator<Double> keyIterator= inputMap.keySet().iterator();
        keyIterator.next();
        double minScale = keyIterator.next() - inputMap.firstKey();

        for (int i = 0; i < curveCount; i++) {
            int offset = i*3;
            if (output == null || output[offset+2] <= minScale || output[offset+1] < -minScale || output[offset] < 0) {
                for (double windowSize = minScale / 10; (output == null || output[offset+2] <= minScale || output[offset+1] < 0) && windowSize <= (minScale / 2); windowSize += minScale / 10) {
                    obs.clear();
                    SortedMap<Double, Double> averaged = MovingAverage.averagedMap(inputMap, windowSize);
                    max = 0;
                    maxLoc = 0;
                    for (Double d : averaged.keySet()) {
                        if (averaged.get(d) > max) {
                            maxLoc = d;
                            max = averaged.get(d);
                        }
                    }
                    if (maxLoc == averaged.firstKey()) {
                        maxLoc = 0.0;
                    }
                    double finalMaxLoc1 = maxLoc;

                    averaged.forEach((key, value) -> {
                        obs.add(key, value);
                        if (key > 2 * finalMaxLoc1) {
                            obs.add(((2 * finalMaxLoc1) - key), value);
                        }
                    });
                    try {
                        output = curveFitter.withCount(curveCount).withMaxIterations(100).fit(obs.toList());
                    } catch (Exception ignored) {}
                }
            }

            if (output == null){
                output = new double[curveCount*3];
                for (int j = 0; j < curveCount; j++) {
                    output[j*3] = 0;
                    output[(j*3)+1] = -1.0;
                    output[(j*3)+2] = inputMap.lastKey();
                }
                return output;
            }

            if (output[offset+2] <= minScale || output[offset+1] < -minScale || output[offset] < 0) {
                output[offset] = 0;
                output[offset+1] = -1.0;
                output[offset+2] = inputMap.lastKey();
            }
        }
        return output;
    }

    private double areaUnderCurve(Gaussian gaussian, SortedMap<Double, Double> map, double mean, double sigma) {

        double auc = 0;
        for (Double d : map.keySet()) {
            if ((mean - (3 * sigma)) < d && d < (mean + (3 * sigma))) {
                auc += gaussian.value(d);
            }
        }
        return auc;
    }

    private double areaUnderCurve(SortedMap<Double, Double> map, double mean, double sigma) {

        double auc = 0;

        for (Double d : map.keySet()) {
            if ((mean - (3 * sigma)) < d && d < (mean + (3 * sigma))) {
                auc += map.get(d);
            }
        }
        return auc;
    }

    private double calcRsquared(){
        final double[] residualsSum = {0};
        final double[] totalSum = {0};
        final double rangeMean = sCorrMap.subMap(gaussFitParameters[1] - (3 * gaussFitParameters[2]), gaussFitParameters[1] + (3 * gaussFitParameters[2])).values().stream().mapToDouble(Double::doubleValue).average().getAsDouble();

        sCorrMap.subMap(gaussFitParameters[1] - (3 * gaussFitParameters[2]), gaussFitParameters[1] + (3 * gaussFitParameters[2])).forEach((key, value) -> {
            residualsSum[0] += Math.pow(value - gaussians.value(key), 2);
            totalSum[0] += Math.pow(value - rangeMean, 2);
        });
        return (1-(residualsSum[0]/totalSum[0]));
    }


}
