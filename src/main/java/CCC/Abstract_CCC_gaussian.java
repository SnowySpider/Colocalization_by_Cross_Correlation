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
package CCC;

import net.imagej.Dataset;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import org.apache.commons.io.FileUtils;
import org.scijava.ItemIO;
import org.scijava.plugin.Parameter;
import org.scijava.table.Table;
import org.scijava.table.Tables;
import utils.CrossCorrelationFunctions;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.*;

/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

public abstract class Abstract_CCC_gaussian <R extends RealType<R>, F extends FloatType> extends Abstract_CCC_base {

    @Parameter(label = "Number of Gaussians to fit:", description = "Values > 1 fit a multi-term sum of Gaussians curve to the data", required=false)
    protected int numGaussians2Fit;
    @Parameter(label = "Generate contribution images?", description = "Generates images that highlight the signal from Image 1 and Image 2 that contributed to the result. Uncheck to use less memory.", required = false)
    protected boolean generateContributionImages;

    @Parameter(type = ItemIO.OUTPUT)
    protected Dataset ContributionOf1, ContributionOf2;

    @Parameter(type = ItemIO.OUTPUT, label = "Correlogram data")
    protected Table<org.scijava.table.Column<Double>, Double> resultsTable;

    //region Time-lapse specific variables for Gaussian fit
    @Parameter (type = ItemIO.OUTPUT, label = "Correlation over time")
    protected Table fullTimeCorrelationTableOut;

    protected ArrayList<LinkedHashMap<String,Double>> fullTimeCorrelationTable;
    //endregion

    protected CrossCorrelationFunctions ccFunctions;

    @Override
    protected void initializePlugin(String[] intermediateNames){
        super.initializePlugin(intermediateNames);
        if(numGaussians2Fit == 0)
            numGaussians2Fit = 1;
        if (generateContributionImages) {
            initializeContributionImages();
            maxStatus += dataset1.getFrames();
        } else if (showIntermediates) {

        }
    }

    protected void initializeContributionImages(){
        ContributionOf1 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset1.getName(), inputAxisTypes);
        ContributionOf1.setAxes(inputCalibratedAxes);
        ContributionOf2 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset2.getName(), inputAxisTypes);
        ContributionOf2.setAxes(inputCalibratedAxes);
    }

    protected void initializeData(RandomAccessibleInterval<FloatType> img1, RandomAccessibleInterval<FloatType> img2, RandomAccessibleInterval<R> mask, double [] inputScale, ImgFactory<R> imgFactory){
        ccFunctions = new CrossCorrelationFunctions(img1, img2, mask, scale, imgFactory);
    }

    protected void fitGaussianCurves(){
        statusService.showStatus(currentStatus++, maxStatus,statusBase + "Fitting gaussian to data");
        try{radialProfiler.correlationData.fitGaussianCurve();}
        catch (NullPointerException e){
            logService.warn("Failed to fit gaussian curve to cross correlation of " + dataset1.getName() + " and " + dataset2.getName() + ", suggesting no correlation between the images.\nAcquired data and intermediate correlation images (if the option was selected) will still be shown. Statistical measures will be set to error values (-1).");
            if(generateContributionImages)
                ++currentStatus;
        }
    }

    protected void generateContributionImages(RandomAccessibleInterval <FloatType> img1, RandomAccessibleInterval <FloatType> img2, Img<FloatType> subtracted, RandomAccessibleInterval <R> gaussianCorrelogramOutput, final RandomAccessibleInterval <R> contribution1, final RandomAccessibleInterval <R> contribution2){
        statusService.showStatus(currentStatus++, maxStatus,statusBase + "Determining channel contributions");
        //gaussModifiedCorr = imgFactory.create(img1);
        Img<FloatType> gaussModifiedCorr = ops.create().img(img1, new FloatType());

        ccFunctions.generateGaussianModifiedCCImage(subtracted, gaussModifiedCorr, radialProfiler.correlationData);

        if(showIntermediates){
            LoopBuilder.setImages(gaussianCorrelogramOutput, gaussModifiedCorr).multiThreaded().forEachPixel((a, b) -> a.setReal(b.get()));
        }

        ccFunctions.calculateContributionImages(img1, img2, gaussModifiedCorr, contribution1, contribution2);
    }

    @Override
    protected void generateFullCorrelationTable(){
        List<HashMap<String,Double>> correlationTableList = new ArrayList<>();
        SortedMap<Double, Double> keyMap;

        keyMap = radialProfiler.correlationData.oCorrMap != null ? radialProfiler.correlationData.oCorrMap : radialProfiler.correlationData.sCorrMap;

        keyMap.keySet().stream().forEachOrdered((d) -> {
            LinkedHashMap<String, Double> row = new LinkedHashMap<String, Double>();
            row.put("Distance (" + getUnitType() +")", (getSigDigits(d)));
            if(radialProfiler.correlationData.oCorrMap != null) row.put("Original CC", getSigDigits(radialProfiler.correlationData.oCorrMap.get(d)));
            if(radialProfiler.correlationData.sCorrMap != null) row.put("Subtracted CC", getSigDigits(radialProfiler.correlationData.sCorrMap.get(d)));
            if(radialProfiler.correlationData.gaussians != null) {
                for (int i = 0; i < numGaussians2Fit; i++) {
                    row.put(("Gaussian fit-" + (i+1)), getSigDigits(radialProfiler.correlationData.gaussians.getGaussian(i).value(d)));
                }
            }
            correlationTableList.add(row);
        });
        correlationTable = Tables.wrap(correlationTableList, null);
    }

    protected void generateResults(){
        //For non-time lapse data, this just returns the Hashmap of the results

        if(dataset1.getFrames() == 1){
            LinkedHashMap<String, Double> gaussResultsHash = new LinkedHashMap<>();
            for (int i = 0; i < numGaussians2Fit; i++) {
                gaussResultsHash.put("Mean"+(i+1)+" (" + getUnitType() + ")", getSigDigits(radialProfiler.correlationData.getGaussianMean(i)));
                gaussResultsHash.put("StDev"+(i+1)+" (" + getUnitType() + ")", getSigDigits(radialProfiler.correlationData.getGaussianSD(i)));
                gaussResultsHash.put("Height"+(i+1), getSigDigits(radialProfiler.correlationData.getGaussianPeakHeight(i)));
                if (radialProfiler.correlationData.hasConfidence()) gaussResultsHash.put(("Confidence"+(i+1)), getSigDigits(radialProfiler.correlationData.getConfidence(i)));
            }
            gaussResultsHash.put("R-squared", getSigDigits(radialProfiler.correlationData.getRSquared()));

            generateResultsTable(gaussResultsHash);
            addGaussianToSummaryFile(gaussResultsHash);
        }
        else{
            LinkedHashMap<String, Double> bestFrame = getBestFrameResults();
            for (int i = 0; i < radialProfiler.correlationData.curveCount; i++) {
                generateResultsTable(bestFrame);
                addGaussianToSummaryFile(bestFrame);
            }
            fullTimeCorrelationTableOut = Tables.wrap(fullTimeCorrelationTable, null);
        }
    }

    protected void generateResultsTable(LinkedHashMap<String, Double> resultsData){
        resultsData.forEach((key, value) ->{value = getSigDigits(value);});
        resultsTable = Tables.wrap(resultsData, null);
    }

    protected void addToFullTimeResultsTable(long frame) {
        if (fullTimeCorrelationTable == null) fullTimeCorrelationTable = new ArrayList<LinkedHashMap<String,Double>>();

        LinkedHashMap<String, Double> gaussianMap = new LinkedHashMap<>();
        for (int i = 0; i < numGaussians2Fit; i++) {
            gaussianMap.put(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "Time (" + calibratedTime.get().unit() + ")" : "Frame",getSigDigits(calibratedTime.get().calibratedValue(frame)));
            gaussianMap.put("Mean-" + (i+1), getSigDigits(radialProfiler.correlationData.getGaussianMean(i)));
            gaussianMap.put("SD-"+ (i+1), getSigDigits(radialProfiler.correlationData.getGaussianSD(i)));
            gaussianMap.put("Gaussian height-"+ (i+1), getSigDigits(radialProfiler.correlationData.getGaussianPeakHeight(i)));
            if (radialProfiler.correlationData.hasConfidence())
                gaussianMap.put("Confidence-"+ (i+1), getSigDigits(radialProfiler.correlationData.getConfidence(i)));
        }
        gaussianMap.put("R-squared", getSigDigits(radialProfiler.correlationData.getRSquared()));

        fullTimeCorrelationTable.add(gaussianMap);
    }

    protected LinkedHashMap<String, Double> getFrameResults(long frame){
        return fullTimeCorrelationTable.get((int)frame);
    }

    protected LinkedHashMap<String, Double> getBestFrameResults(){
        Comparator<LinkedHashMap<String,Double>> tableCompare = (o1, o2) -> {
            Double o1Max = 0.0;
            Double o2Max = 0.0;
            if (o1.get("Confidence-1") != null) {
                for (int i = 0; i < numGaussians2Fit; i++) {
                    o1Max = Math.max(o1Max, o1.get("Confidence-"+(i+1)));
                    o2Max = Math.max(o2Max, o2.get("Confidence-"+(i+1)));
                }
                return o1Max.compareTo(o2Max);
            }
            else{
                for (int i = 0; i < numGaussians2Fit; i++) {
                    o1Max = Math.max(o1Max, o1.get("SD-"+(i+1)));
                    o2Max = Math.max(o2Max, o2.get("SD-"+(i+1)));
                }
                return o2Max.compareTo(o1Max);
            }
        };
        return fullTimeCorrelationTable.stream().max(tableCompare).isPresent() ? fullTimeCorrelationTable.stream().max(tableCompare).get() : null;
    }

    protected void addGaussianToSummaryFile(LinkedHashMap<String, Double> resultsData){
        summary = summary + "Fit a gaussian curve to the cross-correlation of: \n\""
                + dataset1.getName() +
                "\"\n with \n\"" +
                dataset2.getName() +
                "\"\n using the mask \n\"" +
                (maskAbsent? "No mask selected" : maskDataset.getName()) + "\":\n\n";


        resultsData.forEach((key, value) ->{
            summary = summary + key + ": " + getSigDigits(value) + "\n";
        });
        summary = summary + "\n\nFor time-lapse data, only the frame with the highest confidence is shown here. See the Correlation over time data table for complete results.";
        summary = summary + "\n\nConfidence values below 0.1 are considered low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nSee the website for more details.";
        summary = summary + "\n\nAny negative confidence and R-squared values are an error result that indicate that a Gaussian curve could not be fit to the data.\nThis could be because no spatial correlation exists, or because the mask was too narrowly defined. The mask should NOT be a simple segmentation of the objects you want to measure.\nSee the website for details on how to generate a proper mask for your data.";

        if(resultsData.get("R-squared") != null && resultsData.get("R-squared") < 0.05){
            summary = summary + "\n\nThe R-squared value for the gaussian regression is very low.\nThis can indicate a low signal to noise ratio, or that no spatial correlation exists and the curve was fit to image noise.\nSee website for more details.";
        }
    }

    protected void saveResultsToFolder(){
        super.saveResultsToFolder();
        try {
            if(resultsTable != null) ioService.save(resultsTable,saveFolder.getAbsolutePath() + File.separator + "CC Results.csv" );
            FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + File.separator + "Summary.txt"), summary, (Charset) null);
            if(fullTimeCorrelationTable != null){
                ioService.save(fullTimeCorrelationTableOut, saveFolder.getAbsolutePath() + File.separator + "Gaussian fits over time.csv");
            }
            if(generateContributionImages){
                saveDatasetsToFolder(ContributionOf1, ContributionOf2);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

