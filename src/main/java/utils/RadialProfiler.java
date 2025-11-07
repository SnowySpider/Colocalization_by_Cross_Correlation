/*-
 * #%L
 * Scijava plugin for spatial correlation
 * %%
 * Copyright (C) 2019 - 2024 Andrew McCall, University at Buffalo
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

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

import java.util.*;


public class RadialProfiler {

    public CorrelationData correlationData;

    private long[] dimensions;

    private int nDims;

    private double[] scale;

    public RadialProfiler(RandomAccessibleInterval input, double[] inputScale) throws Exception {
        this.initializeToImageDimensions(input, inputScale);
        correlationData = new CorrelationData();
    }

    public RadialProfiler(RandomAccessibleInterval input, double[] inputScale, int curveFitCount) throws Exception {
        this.initializeToImageDimensions(input, inputScale);
        correlationData = new CorrelationData(curveFitCount);
    }

    //This will need to throw an exception, in case the scale an input dimensions don't match
    public void initializeToImageDimensions(RandomAccessibleInterval input, double[] inputScale) throws Exception {
        //get image dimensions and center
        nDims = input.numDimensions();
        if (nDims != inputScale.length)
            throw new Exception("Input image number of dimensions do not match scale number of dimensions");
        scale = inputScale.clone();
        dimensions = new long[nDims];
        input.dimensions(dimensions);

        //BDscale = BigDecimal.valueOf(inputScale[0]).scale() + 3;
    }

    public void calculateOCorrProfile(RandomAccessibleInterval origCorrelation) {
        correlationData.oCorrMap = new TreeMap<>();

        calculateSingleProfile(origCorrelation, correlationData.oCorrMap);
    }

    public void calculateSCorrProfile(RandomAccessibleInterval subtractedCorrelation) {
        correlationData.sCorrMap = new TreeMap<>();

        calculateSingleProfile(subtractedCorrelation, correlationData.sCorrMap);
    }

    public void calculateBothProfiles(RandomAccessibleInterval origCorrelation, RandomAccessibleInterval subtractedCorrelation) {
        correlationData.oCorrMap = new TreeMap<>();
        correlationData.sCorrMap = new TreeMap<>();

        calculateSingleProfile(origCorrelation, correlationData.oCorrMap);
        calculateSingleProfile(subtractedCorrelation, correlationData.sCorrMap);
    }

    private <T extends RealType> void calculateSingleProfile(RandomAccessibleInterval<T> input, Map<Double, Double> output) {
        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = (((double) dimensions[i])-1.0) / 2;
        }

        Map<Double, Double[]> tempMap = Collections.synchronizedMap(new HashMap<>());
        //loop through all points, determine distance (scaled) and bin

        Parallelization.runMultiThreaded(() -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List<Interval> chunks = IntervalChunks.chunkInterval(input, numTasks);

            taskExecutor.forEach(chunks, chunk -> {
                Cursor<T> looper = Views.interval(input, chunk).localizingCursor();
                while (looper.hasNext()) {
                    looper.fwd();
                    double LscaledSq = 0;
                    for (int i = 0; i < nDims; ++i) {
                        LscaledSq += Math.pow((looper.getDoublePosition(i) - center[i]) * scale[i], 2);
                    }
                    Double distance = Math.sqrt(LscaledSq);
                    synchronized (tempMap) {
                        if (tempMap.containsKey(distance)) {
                            tempMap.get(distance)[0] += looper.get().getRealDouble();
                            tempMap.get(distance)[1] += 1;
                        } else {
                            tempMap.put(distance, new Double[2]);
                            tempMap.get(distance)[0] = looper.get().getRealDouble();
                            tempMap.get(distance)[1] = 1.0;
                        }
                    }
                }
            });
        });

        tempMap.forEach((key,value) -> {
            output.put(key.doubleValue(), (value[0]/value[1]));
        });
    }
}
