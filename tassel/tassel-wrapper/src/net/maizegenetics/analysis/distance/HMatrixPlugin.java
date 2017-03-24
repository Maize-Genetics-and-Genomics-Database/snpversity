/*
 *  HMatrixPlugin
 * 
 *  Created on Oct 23, 2015
 */
package net.maizegenetics.analysis.distance;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.taxa.distance.DistanceMatrix;

/**
 * @author Josh Lamos-Sweeney
 * @author Yaw Nti-Addae
 * @author Kelly Robbins
 * @author Terry Casstevens
 */
public class HMatrixPlugin extends AbstractPlugin {

    private PluginParameter<DistanceMatrix> myAMatrix = new PluginParameter.Builder<>("pedigreeMatrix", null, DistanceMatrix.class)
            .description("Pedigree Matrix (A Matrix)")
            .required(true)
            .distanceMatrix()
            .build();

    private PluginParameter<DistanceMatrix> myGMatrix = new PluginParameter.Builder<>("kinshipMatrix", null, DistanceMatrix.class)
            .description("Kinship Matrix (G Matrix)")
            .required(true)
            .distanceMatrix()
            .build();

    private PluginParameter<Double> myWeight = new PluginParameter.Builder<>("weight", 1.0, Double.class)
            .description("Weight")
            .build();

    public HMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> result = new ArrayList<>();

        DoubleMatrix aDoubleMatrix = DoubleMatrixFactory.DEFAULT.make(aMatrix());
        aDoubleMatrix.invert();
        DistanceMatrix aInverse = new DistanceMatrix(aDoubleMatrix.toArray(), aMatrix().getTaxaList());

        DoubleMatrix gDoubleMatrix = DoubleMatrixFactory.DEFAULT.make(gMatrix());
        gDoubleMatrix.invert();
        DistanceMatrix gInverse = new DistanceMatrix(gDoubleMatrix.toArray(), gMatrix().getTaxaList());

        DistanceMatrix matrix = generateCombinedMatrix(aInverse, gInverse, weight());
        result.add(new Datum("Combined A and G Matrix", matrix, null));
        result.add(new Datum("A Matrix Inverse", aInverse, null));
        result.add(new Datum("G Matrix Inverse", gInverse, null));

        return new DataSet(result, this);

    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(HMatrixPlugin.class);
    // }
    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    //public <Type> runPlugin(DataSet input) {
    //    return (<Type>) performFunction(input).getData(0).getData();
    //}
    /**
     * Pedigree Matrix
     *
     * @return Pedigree Matrix
     */
    public DistanceMatrix aMatrix() {
        return myAMatrix.value();
    }

    /**
     * Set Pedigree Matrix. Pedigree Matrix
     *
     * @param value Pedigree Matrix
     *
     * @return this plugin
     */
    public HMatrixPlugin aMatrix(DistanceMatrix value) {
        myAMatrix = new PluginParameter<>(myAMatrix, value);
        return this;
    }

    /**
     * Kinship Matrix
     *
     * @return Kinship Matrix
     */
    public DistanceMatrix gMatrix() {
        return myGMatrix.value();
    }

    /**
     * Set Kinship Matrix. Kinship Matrix
     *
     * @param value Kinship Matrix
     *
     * @return this plugin
     */
    public HMatrixPlugin gMatrix(DistanceMatrix value) {
        myGMatrix = new PluginParameter<>(myGMatrix, value);
        return this;
    }

    /**
     * Weight
     *
     * @return Weight
     */
    public Double weight() {
        return myWeight.value();
    }

    /**
     * Set Weight. Weight
     *
     * @param value Weight
     *
     * @return this plugin
     */
    public HMatrixPlugin weight(Double value) {
        myWeight = new PluginParameter<>(myWeight, value);
        return this;
    }

    /**
     * Given the inverse of an A matrix (Pedigree Kinship matrix) and a G Matrix
     * (Genetic Kinship Matrix) which contains mostly entries in the A Matrix
     * creates a combined (H) matrix.
     *
     * @param aInverse Inverted Pedigree kinship(A) matrix
     * @param gInverse Inverted kinship(G) matrix. Should contain mostly Taxa in
     * common with A matrix for best results
     * @param weight 0.0 to 1.0, the 'Trust' put in the G matrix. When creating
     * the combined matrix, uses inverse G * weight + inverse A * 1-weight
     * @return Combined matrix based on aInverse and gInverse
     */
    // Shorthand notation used: A' = A^-1 (A inverse) = A prime
    // a[1,1] is a11 from the research paper
    public DistanceMatrix generateCombinedMatrix(DistanceMatrix aInverse, DistanceMatrix gInverse, double weight) {
        List<Taxon> aTaxa = aInverse.getTaxaList();
        List<Taxon> gTaxa = gInverse.getTaxaList();

        //Generate Taxa intersections
        List<Taxon> aOnlyTaxa = new LinkedList<>(aTaxa);
        aOnlyTaxa.removeAll(gTaxa);

        List<Taxon> gOnlyTaxa = new LinkedList<>(gTaxa);
        gOnlyTaxa.removeAll(aTaxa);

        List<Taxon> unionTaxa = new LinkedList<>(aTaxa);
        unionTaxa.removeAll(aOnlyTaxa);

        //Entries not in A are added to the end of A as identity.
        //A[1]
        List<Taxon> matrixOrder = new ArrayList<>(aOnlyTaxa);
        matrixOrder.addAll(gOnlyTaxa);
        //A[2]
        matrixOrder.addAll(unionTaxa);

        double[][] aPrime = aInverse.getDistances();
        double[][] gPrime = gInverse.getDistances();

        int outputSize = aTaxa.size() + gOnlyTaxa.size();
        double[][] hPrime = new double[outputSize][outputSize];

        //Combined Matrix entries for A[1][1] are copied from A matrix
        for (int i = 0; i < aPrime.length; i++) {
            int newI = matrixOrder.indexOf(aTaxa.get(i));
            for (int j = 0; j < aPrime.length; j++) {
                int newJ = matrixOrder.indexOf(aTaxa.get(j));
                hPrime[newI][newJ] = aPrime[i][j];
            }
        }
        //Add Taxa only in the GMatrix to the A Matrix with an identity of 1 and no other values
        for (Taxon gOnly : gOnlyTaxa) {
            int newI = matrixOrder.indexOf(gOnly);
            hPrime[newI][newI] = 1.0;
        }

        //Now that we have a combined matrix of
        //a11 a12
        //a21 a22
        //we can replace a22 with
        //w*g - (1-w)a22
        for (int i = 0; i < gTaxa.size(); i++) {
            int newI = matrixOrder.indexOf(gTaxa.get(i));
            for (int j = 0; j < gTaxa.size(); j++) {
                int newJ = matrixOrder.indexOf(gTaxa.get(j));
                hPrime[newI][newJ] = weight * gPrime[i][j] - (1.0 - weight) * hPrime[newI][newJ];//TODO discuss weighting. 1 becomes naive
            }
        }

        for (int i = aOnlyTaxa.size(); i < aTaxa.size(); i++) {
        }
        TaxaListBuilder builder = new TaxaListBuilder();
        for (Taxon current : matrixOrder) {
            builder.add(current);
        }
        DistanceMatrix hInverse = new DistanceMatrix(hPrime, builder.build());
        return hInverse;
    }

    /**
     * Gets the Eigenvalue Decomposition of a distance matrix. Helpful in
     * determining the optimal weights of a combined matrix.
     *
     * @param input Input matrix
     * @return Eigenvalue Decomposition of the matrix
     */
    private static EigenvalueDecomposition decompose(DistanceMatrix input) {
        return DoubleMatrixFactory.DEFAULT.make(input.getDistances()).getEigenvalueDecomposition();
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = HMatrixPlugin.class.getResource("/net/maizegenetics/analysis/images/hmatrix.png");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Combined A and G Relationship Matrix";
    }

    @Override
    public String getToolTipText() {
        return "Create Combined A and G Relationship Matrix (H Matrix)";
    }

    @Override
    public String getCitation() {
        return "Lamos-Sweeney J, Nti-Addae Y, Robbins K, Casstevens T. (Oct. 2015) Second Tassel Hackathon.";
    }

}
