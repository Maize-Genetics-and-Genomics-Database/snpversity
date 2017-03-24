// DistanceMatrixUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.taxa.distance;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;

import java.io.Serializable;
import java.util.Arrays;


/**
 * Auxillary functions for distance matrices<p>
 *
 * @version $Id: DistanceMatrixUtils.java,v 1.1 2007/01/12 03:26:14 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public class DistanceMatrixUtils implements Serializable {

	/**
	 * compute squared distance to second distance matrix.
	 * If both matrices have the same size it is assumed that the order of the taxa
	 * is identical.
	 */
	public static double squaredDistance(DistanceMatrix mat1, DistanceMatrix mat2, boolean weighted) {

		boolean aliasNeeded = false;
		if (mat1.getSize() != mat2.getSize())
		{
			aliasNeeded = true;
		}

		int[] alias = null;

		if (aliasNeeded) {
			if (mat1.getSize() > mat2.getSize()) {
				//swap so mat1 is the smaller of the two
				DistanceMatrix temp = mat2;
				mat2 = mat1;
				mat1 = temp;
			}
			alias = new int[mat1.getSize()];
			for (int i = 0; i < alias.length; i++) {
				alias[i] = mat2.whichIdNumber(mat1.getTaxon(i).getName());
			}
		} else {
			alias = new int[mat1.getSize()];
			for (int i = 0; i < alias.length; i++) {
				alias[i] = i;
			}
		}


		double sum = 0;
		int ai;
		final double[][] mat1Distance = mat1.getDistances();
		final double[][] mat2Distance = mat2.getDistances();
		for (int i = 0; i < mat1.getSize()-1; i++)
		{
			ai = alias[i];

			for (int j = i+1; j < mat1.getSize(); j++)
			{
				double diff = mat1Distance[i][j] - mat2Distance[ai][alias[j]];
				double weight;
				if (weighted)
				{
					// Fitch-Margoliash weight
					// (variances proportional to distances)
					weight = 1.0/(mat1Distance[i][j]*mat2Distance[ai][alias[j]]);
				}
				else
				{
					// Cavalli-Sforza-Edwards weight
					// (homogeneity of variances)
					weight = 1.0;
				}
				sum += weight*diff*diff;
			}
		}

		return 2.0*sum; // we counted only half the matrix
	}

	/**
	 * Returns a distance matrix with the specified taxa removed.
	 */
	public static DistanceMatrix minus(DistanceMatrix parent, int taxaToRemove) {

		int size = parent.numberOfTaxa() - 1;

		double[][] distances = new double[size][size];
		Taxon[] ids = new Taxon[size];
		int counti = 0, countj = 0;
		for (int i = 0; i < size; i++) {
			if (counti == taxaToRemove) {
				counti += 1;
			}
			ids[i] = parent.getTaxon(counti);

			countj = 0;
			final double[][] parentDistance = parent.getDistances();
			for (int j = 0; j < size; j++) {
				if (countj == taxaToRemove) {
					countj += 1;
				}
				distances[i][j] = parentDistance[counti][countj];
				countj += 1;
			}
			counti += 1;
		}
        TaxaList tl=new TaxaListBuilder().addAll(ids).build();
		DistanceMatrix smaller = new DistanceMatrix(distances, tl);

		return smaller;
	}
        
    /**
     * Calculates the IBS distance between two taxa with bitsets for for major and minor allele
     * @param iMajor
     * @param iMinor
     * @param jMajor
     * @param jMinor
     * @return 
     */ 
    public static double getIBSDistance(long[] iMajor, long[] iMinor, long[] jMajor, long[] jMinor) {
        int sameCnt=0, diffCnt=0, hetCnt=0;
        for(int x=0; x<iMajor.length; x++) {
            long same=(iMajor[x]&jMajor[x])|(iMinor[x]&jMinor[x]);
            long diff=(iMajor[x]&jMinor[x])|(iMinor[x]&jMajor[x]);
            long hets=same&diff;
            sameCnt+=BitUtil.pop(same);
            diffCnt+=BitUtil.pop(diff);
            hetCnt+=BitUtil.pop(hets);
        }
        double identity=(double)(sameCnt+(hetCnt/2))/(double)(sameCnt+diffCnt+hetCnt);
        double dist=1-identity;
        return dist;
    }
    
    public static double getIBSDistance(BitSet iMajor, BitSet iMinor, BitSet jMajor, BitSet jMinor) {
        return getIBSDistance(iMajor.getBits(), iMinor.getBits(), jMajor.getBits(), jMinor.getBits());
    }

    /**
     * @param parent	the DistanceMatrix from which to extract a subset
     * @param taxaToKeep	an index of the taxa to keep
     * @return A DistanceMatrix with all the taxa that are in both parent and taxaToKeep in the same order as taxaToKeep
     */
    public static DistanceMatrix keepTaxa(DistanceMatrix parent, int[] taxaToKeep) {
    	int ntaxa = taxaToKeep.length;
    	double[][] newDistances = new double[ntaxa][ntaxa];
    	for (int r = 0; r < ntaxa; r++) {
    		for (int c = 0; c < ntaxa; c++) {
    			newDistances[r][c] = parent.getDistance(taxaToKeep[r], taxaToKeep[c]);
    		}
    	}

    	TaxaListBuilder taxaBuilder = new TaxaListBuilder();
    	for (int ndx : taxaToKeep) taxaBuilder.add(parent.getTaxon(ndx));
    	TaxaList taxaListToKeep = taxaBuilder.build();

    	return new DistanceMatrix(newDistances, taxaListToKeep);
    }
    
    /**
     * @param parent	the DistanceMatrix from which to extract a subset
     * @param taxaToKeep	a TaxaList of taxa to be kept 
     * @return	a DistanceMatrix that contains only the taxa that are in both taxaToKeep and parent.
     * 	The taxa will be in the same order as taxaToKeep.
     */
    public static DistanceMatrix keepTaxa(DistanceMatrix parent, TaxaList taxaToKeep) {
    	int[] keepIndex = taxaToKeep.stream()
    		.mapToInt(t -> parent.whichIdNumber(t))
    		.filter(i -> i > -1)
    		.toArray();
    	return keepTaxa(parent, keepIndex);
    }
}
