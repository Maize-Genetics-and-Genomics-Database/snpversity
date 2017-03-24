package net.maizegenetics.taxa.distance;


import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.TableReport;

import java.io.Serializable;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Jan 10, 2007
 * Time: 12:13:28 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TaxaListMatrix extends Serializable, TableReport {

    /**
     * Return TaxaList of this Matrix.
     */
    public TaxaList getTaxaList();

    /** returns representation of this alignment as a string */
    String toString();

    /**
	 * Returns the number of rows and columns that the distance matrix has.
	 */
    int getSize();

    /**
	 * Returns the distances as a 2-dimensional array of doubles. Matrix is cloned first so it can be altered freely.
	 */
    double[][] getClonedDistances();

    /*
    Returns the distance calculated for two taxa with the indices row and col
     */
    double getDistance(int row, int col);

    /**
	 * Returns the mean pairwise distance of this matrix
	 */
    double meanDistance();

    	/**
	 * test whether this matrix is a symmetric distance matrix
	 *
	 */
    boolean isSymmetric();
}
