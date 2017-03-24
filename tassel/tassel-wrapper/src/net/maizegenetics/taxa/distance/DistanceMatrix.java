// DistanceMatrix.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.taxa.distance;

import net.maizegenetics.util.FormattedOutput;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * storage for pairwise distance matrices.<p>
 *
 * features: - printing in in PHYLIP format, - computation of (weighted) squared
 * distance to other distance matrix - Fills in all of array...
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class DistanceMatrix implements TaxaListMatrix, TableReport {

    private final TaxaList myTaxaList;
    /**
     * distances [seq1][seq2]
     */
    private double[][] distance = null;
    static final long serialVersionUID = 4725925229860707633L;

    /**
     * constructor taking distances array and IdGroup
     */
    public DistanceMatrix(double[][] distance, TaxaList taxaList) {
        super();
        this.distance = distance;
        myTaxaList = taxaList;
    }

    /**
     * constructor that takes a distance matrix and clones the distances and
     * IdGroup
     */
    public DistanceMatrix(DistanceMatrix dm) {
        double[][] copy = new double[dm.distance.length][];
        for (int i = 0; i < copy.length; i++) {
            copy[i] = new double[dm.distance[i].length];
            System.arraycopy(dm.distance[i], 0, copy[i], 0, dm.distance[i].length);
        }
        distance = copy;
        myTaxaList = dm.myTaxaList;
    }

    /**
     * constructor that takes a distance matrix and clones the distances, of a
     * the identifiers in taxaList.
     */
    public DistanceMatrix(DistanceMatrix dm, TaxaList subset) {

        int index1, index2;

        distance = new double[subset.numberOfTaxa()][subset.numberOfTaxa()];
        for (int i = 0; i < distance.length; i++) {
            index1 = dm.whichIdNumber(subset.taxaName(i));
            distance[i][i] = dm.distance[index1][index1];
            for (int j = 0; j < i; j++) {
                index2 = dm.whichIdNumber(subset.taxaName(j));
                distance[i][j] = dm.distance[index1][index2];
                distance[j][i] = distance[i][j];
            }
        }
        myTaxaList = subset;
    }

    /**
     * print alignment (PHYLIP format)
     */
    public void printPHYLIP(PrintWriter out) throws IOException {
        // PHYLIP header line
        out.println("  " + distance.length);
        FormattedOutput format = FormattedOutput.getInstance();

        for (int i = 0; i < distance.length; i++) {
            format.displayLabel(out,
                    myTaxaList.taxaName(i), 10);
            out.print("      ");

            for (int j = 0; j < distance.length; j++) {
                // Chunks of 6 blocks each
                if (j % 6 == 0 && j != 0) {
                    out.println();
                    out.print("                ");
                }

                out.print("  ");
                format.displayDecimal(out, distance[i][j], 5);
            }
            out.println();
        }
    }

    /**
     * returns representation of this alignment as a string
     */
    @Override
    public String toString() {

        StringWriter sw = new StringWriter();
        try {
            printPHYLIP(new PrintWriter(sw));
        } catch (Exception e) {
            e.printStackTrace();
        }

        return sw.toString();
    }

    /**
     * compute squared distance to second distance matrix
     */
    public double squaredDistance(DistanceMatrix mat, boolean weighted) {
        double sum = 0;
        for (int i = 0; i < distance.length - 1; i++) {
            for (int j = i + 1; j < distance.length; j++) {
                double diff = distance[i][j] - mat.distance[i][j];
                double weight;
                if (weighted) {
                    // Fitch-Margoliash weight
                    // (variances proportional to distances)
                    weight = 1.0 / (distance[i][j] * distance[i][j]);
                } else {
                    // Cavalli-Sforza-Edwards weight
                    // (homogeneity of variances)
                    weight = 1.0;
                }
                sum += weight * diff * diff;
            }
        }

        return 2.0 * sum; // we counted only half the matrix
    }

    /**
     * compute absolute distance to second distance matrix
     */
    public double absoluteDistance(DistanceMatrix mat) {
        double sum = 0;
        for (int i = 0; i < distance.length - 1; i++) {
            for (int j = i + 1; j < distance.length; j++) {
                double diff
                        = Math.abs(distance[i][j] - mat.distance[i][j]);

                sum += diff;
            }
        }

        return 2.0 * sum; // we counted only half the matrix
    }

    /**
     * Returns the number of rows and columns that the distance matrix has.
     */
    public int getSize() {
        return distance.length;
    }

    /**
     * Returns the distances as a 2-dimensional array of doubles. Matrix is
     * cloned first so it can be altered freely.
     */
    public final double[][] getClonedDistances() {
        double[][] copy = new double[distance.length][];
        for (int i = 0; i < copy.length; i++) {
            copy[i] = new double[distance[i].length];
            System.arraycopy(distance[i], 0, copy[i], 0, distance[i].length);
        }
        return copy;
    }

    /**
     * Returns the distances as a 2-dimensional array of doubles (in the actual
     * array used to store the distances)
     */
    public final double[][] getDistances() {
        return getClonedDistances();
    }

    public final double getDistance(final int row, final int col) {
        return distance[row][col];
    }

    /**
     * Sets both upper and lower triangles.
     *
     * @deprecated Needs to have a Builder
     */
    public void setDistance(int i, int j, double dist) {
        distance[i][j] = distance[j][i] = dist;
    }

    /**
     * Returns the mean pairwise distance of this matrix
     */
    @Override
    public double meanDistance() {
        double dist = 0.0;
        int count = 0;
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                if ((i != j) && (!Double.isNaN(distance[i][j]))) {
                    dist += distance[i][j];
                    count += 1;
                }
            }
        }
        return dist / (double) count;
    }

    //IdGroup interface
    public Taxon getTaxon(int i) {
        return myTaxaList.get(i);
    }

    public int numberOfTaxa() {
        return myTaxaList.numberOfTaxa();
    }

    public int whichIdNumber(String name) {
        return myTaxaList.indexOf(name);
    }

    public int whichIdNumber(Taxon id) {
        return myTaxaList.indexOf(id);
    }

    /**
     * Return TaxaList of this alignment.
     */
    public TaxaList getTaxaList() {
        return myTaxaList;
    }

    /**
     * test whether this matrix is a symmetric distance matrix
     *
     */
    @Override
    public boolean isSymmetric() {
        for (int i = 0; i < distance.length; i++) {
            if (distance[i][i] != 0) {
                return false;
            }
        }
        for (int i = 0; i < distance.length - 1; i++) {
            for (int j = i + 1; j < distance.length; j++) {
                if (distance[i][j] != distance[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }

    private boolean isIn(int value, int[] set) {
        if (set == null) {
            return false;
        }
        for (int i = 0; i < set.length; i++) {
            if (set[i] == value) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param fromIndex the index of the thing (taxa,sequence) from which we
     * want to find the closest (excluding self)
     * @param exclusion indexes of things that should not be considered, may be
     * null
     * @return the index of the member closes to the specified
     */
    public int getClosestIndex(int fromIndex, int[] exclusion) {
        double min = Double.POSITIVE_INFINITY;
        int index = -1;
        for (int i = 0; i < distance.length; i++) {
            if (i != fromIndex && !isIn(i, exclusion)) {
                double d = distance[fromIndex][i];
                if (d < min) {
                    min = d;
                    index = i;
                }
            }
        }
        return index;
    }

    public static DistanceMatrix hadamardProduct(DistanceMatrix m0, DistanceMatrix m1) {
        int n = m0.distance.length;
        if (m0.distance.length != n) 
            throw new IllegalArgumentException("Matrices must be of the same dimensions to compute a Hadamard product.");
        
        double[][] product = new double[n][n];
        for (int r = 0; r < n; r++) {
            product[r][r] = m0.distance[r][r] * m0.distance[r][r];
            for (int c = r + 1; c < n; c++) {
                product[r][c] = product[c][r] = m0.distance[r][c] * m0.distance[r][c];
            }
        }
        
        return new DistanceMatrix(product, m0.getTaxaList());    
    }
    
    @Override
    public Object[] getTableColumnNames() {
        String[] colNames = new String[getSize() + 1];
        colNames[0] = "Taxa";
        for (int i = 0; i < distance[0].length; i++) {
            colNames[i + 1] = getTaxon(i).toString();
        }
        return colNames;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    @Override
    public Object[] getRow(long rowLong) {

        int row = (int) rowLong;
        Object[] result = new Object[distance[row].length + 1];
        result[0] = getTaxon(row);
        for (int j = 1; j <= distance[row].length; j++) {
            result[j] = "" + distance[row][j - 1];
        }

        return result;

    }

    @Override
    public String getTableTitle() {
        return "Distance Matrix";
    }

    @Override
    public long getRowCount() {
        if (distance != null) {
            return distance.length;
        } else {
            return 0;
        }
    }

    @Override
    public long getElementCount() {
        return getRowCount() * getColumnCount();
    }

    @Override
    public int getColumnCount() {
        if ((distance != null) && distance[0] != null) {
            return distance[0].length + 1;
        } else {
            return 0;
        }
    }

    @Override
    public Object getValueAt(long rowIndex, int columnIndex) {
        if (columnIndex == 0) {
            return getTaxon((int) rowIndex);
        }
        return distance[(int) rowIndex][columnIndex - 1];
    }

    public String getColumnName(int col) {
        if (col == 0) {
            return "Taxa";
        }
        return getTaxon(col - 1).toString();
    }

}
