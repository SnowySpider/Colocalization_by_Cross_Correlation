/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package utils;

import java.util.Arrays;

import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.DifferentiableUnivariateFunction;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

/**
 * <a href="http://en.wikipedia.org/wiki/Gaussian_function">
 *  Gaussian</a> function.
 *
 * @since 3.0
 */
public class nGaussian implements UnivariateDifferentiableFunction, DifferentiableUnivariateFunction {

    private final int count;

    private final Gaussian[] gaussians;

    /**
     * Gaussian with given normalization factor, mean and standard deviation.
     *
     * @param gaussians Abstract array of org.apache.commons.math3.analysis.function.Gaussian classes that make up this sum
     */
    public nGaussian(Gaussian ... gaussians) {
        this.count = gaussians.length;
        this.gaussians = gaussians.clone();
    }

    public nGaussian(double ... gaussianParams){
        if(gaussianParams.length%3 != 0){
            throw new DimensionMismatchException(gaussianParams.length, gaussianParams.length + gaussianParams.length%3);
        }
        this.count = gaussianParams.length/3;
        gaussians = new Gaussian[this.count];
        int i = 0;
        for (int j = 0; j < this.count; j++) {
            gaussians[j] = new Gaussian(gaussianParams[i++],gaussianParams[i++],gaussianParams[i++]);
        }
    }

    public void setGaussian(int index, Gaussian inputGaussian){
        gaussians[index] = inputGaussian;
    }

    public Gaussian getGaussian(int index){
        return gaussians[index];
    }

    /** {@inheritDoc} */
    public double value(double x) {
        double returnValue = 0.0;
        for(Gaussian gaussian: gaussians)
            returnValue += gaussian.value(x);
        return returnValue;
    }

    //Does not work ; was not updated. Simply needed for abstract implementation
    @Deprecated
    public UnivariateFunction derivative() {
        return FunctionUtils.toDifferentiableUnivariateFunction(this).derivative();
    }

    /**
     * Parametric function where the input array contains the parameters of
     * the Gaussian, ordered as follows:
     * <ul>
     *  <li>Norm</li>
     *  <li>Mean</li>
     *  <li>Standard deviation</li>
     * </ul>
     */
    public static class Parametric implements ParametricUnivariateFunction {
        /**
         * Computes the value of the sum of n-Gaussians at {@code x}.
         *
         * @param x Value for which the function must be computed.
         * @param param Values of norm, mean and standard deviation for each Gaussian.
         * @return the value of the function.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws DimensionMismatchException if the size of {@code param} is
         * not 3.
         * @throws NotStrictlyPositiveException if {@code param[2]} is negative.
         */
        public double value(double x, double ... param)
                throws NullArgumentException,
                DimensionMismatchException,
                NotStrictlyPositiveException {
            return nGaussian.value(x, param);
        }

        /**
         * Computes the value of the gradient at {@code x}.
         * The components of the gradient vector are the partial
         * derivatives of the function with respect to each of the
         * <em>parameters</em> (norm, mean and standard deviation).
         *
         * @param x Value at which the gradient must be computed.
         * @param param Values of norm, mean and standard deviation.
         * @return the gradient vector at {@code x}.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws DimensionMismatchException if the size of {@code param} is
         * not 3.
         * @throws NotStrictlyPositiveException if {@code param[2]} is negative.
         */
        public double[] gradient(double x, double ... param)
                throws NullArgumentException,
                DimensionMismatchException,
                NotStrictlyPositiveException {
            validateParameters(param);

            //returned double[] should match incoming param number! yay!

            double [] gradients = new double[param.length];

            //This treats each gaussian separately, not sure if that's correct
            for (int i = 0; i < param.length; i += 3) {
                final double norm = param[i];
                final double diff = x - param[i+1];
                final double sigma = param[i+2];
                final double i2s2 = 1 / (2 * sigma * sigma);

                gradients[i] = nGaussian.value(diff, 1, i2s2); //n
                gradients[i+1] = norm * gradients[i] * 2 * i2s2 * diff; //m
                gradients[i+2] = gradients[i+1] * diff / sigma; //s
            }

            return gradients;
        }

        /**
         * Validates parameters to ensure they are appropriate for the evaluation of
         * the {@link #value(double,double[])} and {@link #gradient(double,double[])}
         * methods.
         *
         * @param param Values of norm, mean and standard deviation.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws DimensionMismatchException if the size of {@code param} is
         * not 3.
         * @throws NotStrictlyPositiveException if {@code param[2]} is negative.
         */

    }

    private static double value(double x, double ... param){
        validateParameters(param);
        double sumV = 0;
        int nGauss = (param.length/3);

        double [][] pass = new double[nGauss][3];

        int pAccessor = 0;
        for (int i = 0; i < nGauss; ++i) {
            pass[i][0] = param[pAccessor++];
            pass[i][1] = param[pAccessor++];
            pass[i][2] = param[pAccessor++];
        }
        try {

            for (int i = 0; i < nGauss; i++) {
                final double diff = x - pass[i][1];
                final double i2s2 = 1 / (2 * pass[i][2] * pass[i][2]);
                sumV += nGaussian.singleValue(diff, pass[i][0], i2s2);
            }
        } catch (NotStrictlyPositiveException e) { // NOPMD
            // Do nothing.
        }
        return sumV;
    }
    private static double singleValue(double xMinusMean,
                                double norm,
                                double i2s2) {
        return norm * FastMath.exp(-xMinusMean * xMinusMean * i2s2);
    }

    private static void validateParameters(double[] param)
            throws NullArgumentException,
            DimensionMismatchException,
            NotStrictlyPositiveException {
        if (param == null) {
            throw new NullArgumentException();
        }
//        if (param.length%3!=0) {
//            throw new DimensionMismatchException(param.length, 3);
//        }
        for (int i = 2; i < param.length; i += 3) {
            if (param[i] <= 0) {
                throw new NotStrictlyPositiveException(param[i]);
            }
        }

    }

    /** {@inheritDoc}
     * @since 3.1
     */
    public DerivativeStructure value(final DerivativeStructure t)
            throws DimensionMismatchException {

        try {
            for(Gaussian gaussian: gaussians){
                t.add(gaussian.value(t));
            }
        }
        catch (DimensionMismatchException e) {
            throw e;
        }
        return t;
    }

}