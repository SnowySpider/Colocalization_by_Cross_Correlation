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
package Contribution;

import io.scif.config.SCIFIOConfig;
import io.scif.services.DatasetIOService;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.CalibratedAxis;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetOne;
import net.imglib2.view.Views;
import org.scijava.ItemIO;
import org.scijava.ItemVisibility;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.ContextCommand;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import utils.Contributions;
import utils.CorrelationData;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

@Plugin(type = Command.class, headless = true, menuPath = "Analyze>Colocalization>Colocalization by Cross Correlation>Custom contribution images")
public class Custom_Contribution extends ContextCommand {

    @Parameter
    protected LogService logService;

    @Parameter
    protected StatusService statusService;

    @Parameter
    protected DatasetService datasetService;

    @Parameter
    protected DatasetIOService datasetIOService;

    @Parameter
    protected OpService ops;

    @Parameter(label = "Image 1: ", description = "Image order must persist if using cross-correlation image", persist = false)
    protected Dataset dataset1;

    @Parameter(label = "Image 2: ", persist = false)
    protected Dataset dataset2;

    @Parameter(label = "Use cross-correlation image?", description = "When checked, uses provided cross-correlation image to generate the contribution images. Otherwise generates a normalized Gaussian curve.", required = false)
    protected boolean useCCimage;

    @Parameter(label = "Cross-correlation image: ", description = "The cross-correlation image used to generate the contribution images. Recommend using the 'subtracted cross-correlation' result from CCC", required = false, persist = false)
    protected Dataset ccImage;

    @Parameter(label="Parameters must be given as a comma separated list in order of Height, Mean, SD. Two Gaussian Example:\n"+
            " 3.205E7, 0.7433, 0.3369, 4.841E7, 2.046, 1.342", visibility = ItemVisibility.MESSAGE, required = false)
    protected String msg;

    @Parameter(label="Gaussian curve(s) parameters: ", description = "Comma separated list of Gaussian parameters in repeating Height, Mean, SD order", validater = "checkGaussians")
    protected String values;

    @Parameter(label = "Output directory (leave blank for none):", description = "The directory to automatically save all generated output, including the intermediate images if the \"Show Intermediates\" box is checked", required = false, style="directory")
    protected File saveFolder;

    @Parameter(type = ItemIO.OUTPUT)
    protected Dataset ContributionOf1, ContributionOf2;

    CorrelationData correlationData;

    double [] scale;
    double [] gaussParams;

    protected long[] minDims;
    protected long[] maxDims;
    protected int timeAxis;
    protected SCIFIOConfig config;
    protected CalibratedAxis[] inputCalibratedAxes;
    protected AxisType[] inputAxisTypes;

    @SuppressWarnings("unused")
    public void checkGaussians(){
        this.gaussParams = Arrays.stream(values.split(",")).mapToDouble(Double::parseDouble).toArray();
        if(gaussParams.length % 3 != 0){
            cancel("Gaussian parameters must be provided in triplets.");
        }
        for (int i = 2; i < gaussParams.length; i +=3) {
            if(gaussParams[i] <= 0){
                cancel("All standard deviation values must be > 0.");
            }
        }
    }

    @Override
    public boolean isCanceled() {
        return super.isCanceled();
    }

    @Override
        public void run() {

            if(isCanceled()){
                logService.error(getCancelReason());
                return;
            }

            //region Setup
            if(!useCCimage){
                    ccImage = dataset1.duplicateBlank();
                    LoopBuilder.setImages(ccImage).multiThreaded().forEachPixel(SetOne::setOne);
            }

            correlationData = new CorrelationData(gaussParams);

            if(saveFolder != null) saveFolder.mkdirs();

            config = new SCIFIOConfig();
            config.writerSetFailIfOverwriting(false);

            inputCalibratedAxes = new CalibratedAxis[dataset1.numDimensions()];
            inputAxisTypes = new AxisType[dataset1.numDimensions()];
            for (int i = 0; i < dataset1.numDimensions(); ++i) {
                    inputCalibratedAxes[i] = dataset1.axis(i);
                    inputAxisTypes[i] = dataset1.axis(i).type();
            }

            if(dataset1.getFrames() == 1){
                    scale = new double[dataset1.numDimensions()];
                    for (int i = 0; i < scale.length; i++) {
                            scale[i] = dataset1.averageScale(i);
                    }
            }
            else {
                    timeAxis = dataset1.dimensionIndex(Axes.TIME);
                    minDims = new long[dataset1.numDimensions()];
                    maxDims = dataset1.dimensionsAsLongArray();

                    scale = new double[dataset1.numDimensions()-1];
                    int j = 0;
                    for (int i = 0; i < dataset1.numDimensions(); i++) {
                            if(i != timeAxis) {
                                    maxDims[i] = maxDims[i]-1;
                                    scale[j++] = dataset1.averageScale(i);
                            }
                    }
            }

            ContributionOf1 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset1.getName(), inputAxisTypes);
            ContributionOf1.setAxes(inputCalibratedAxes);
            ContributionOf2 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset2.getName(), inputAxisTypes);
            ContributionOf2.setAxes(inputCalibratedAxes);
            //endregion

            //region single-frame calculation
            if(dataset1.getFrames() == 1){
                    Img<FloatType> gaussModifiedCorr = ops.create().img(dataset1, new FloatType());
                    Contributions.generateGaussianModifiedCCImage(ccImage, gaussModifiedCorr, correlationData, scale);
                    Contributions.calculateContributionImages(dataset1, dataset2, gaussModifiedCorr, ContributionOf1, ContributionOf2, dataset1.getImgPlus().factory());
            }
            //endregion
            //region multi-frame calculation
            else{
                    for (long i = 0; i < dataset1.getFrames(); i++) {
                            setActiveFrame(i);
                            RandomAccessibleInterval temp1 = getActiveFrame(dataset1);
                            RandomAccessibleInterval temp2 = getActiveFrame(dataset2);
                            RandomAccessibleInterval ccImageTemp = getActiveFrame(ccImage);
                            Img<FloatType> gaussModifiedCorr = ops.create().img(ccImageTemp, new FloatType());

                            //Have to do this for every frame for
                            Contributions.generateGaussianModifiedCCImage(ccImageTemp, gaussModifiedCorr, correlationData, scale);

                            Contributions.calculateContributionImages(temp1, temp2, gaussModifiedCorr, getActiveFrame(ContributionOf1), getActiveFrame(ContributionOf2), dataset1.getImgPlus().factory());
                    }
            }
            //endregion

            if(saveFolder != null && !saveFolder.getPath().equals("")) {
                    try {
                            datasetIOService.save(ContributionOf1, saveFolder.getAbsolutePath() + File.separator + ContributionOf1.getName() + ".tif", config);
                            datasetIOService.save(ContributionOf2, saveFolder.getAbsolutePath() + File.separator + ContributionOf2.getName() + ".tif", config);
                    } catch (IOException e) {
                            e.printStackTrace();
                    }
            }

        }

        protected void setActiveFrame(long frame){
                statusService.showProgress((int)frame, (int)dataset1.getFrames());

                minDims[timeAxis] = frame;
                maxDims[timeAxis] = frame;
        }

        protected RandomAccessibleInterval<RealType<?>> getActiveFrame(Dataset in){
                return Views.dropSingletonDimensions(Views.interval(in, minDims, maxDims));
        }
}
