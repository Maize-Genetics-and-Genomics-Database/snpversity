/*
 *  Byte2DBuilder
 */
package net.maizegenetics.dna.snp.byte2d;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

/**
 * Builder to store 2 dimensional byte encoded values.
 *
 * @author Terry Casstevens
 */
public class Byte2DBuilder {

    private SuperByteMatrix myValues = null;
    private final boolean myIsHDF5;
    private IHDF5Writer myHDF5Writer = null;
    private final int myNumSites;
    private int myNumTaxa = 0;
    private final SiteScore.SITE_SCORE_TYPE mySiteScoreType;
    private final TaxaList myTaxaList;

    private Byte2DBuilder(IHDF5Writer writer, int numSites, SiteScore.SITE_SCORE_TYPE siteScoreType, TaxaList taxaList) {
        myIsHDF5 = true;
        myHDF5Writer = writer;
        myNumSites = numSites;
        mySiteScoreType = siteScoreType;
        myTaxaList = taxaList;
    }

    private Byte2DBuilder(int numTaxa, int numSites, SiteScore.SITE_SCORE_TYPE siteScoreType, TaxaList taxaList) {
        myIsHDF5 = false;
        myHDF5Writer = null;
        myNumSites = numSites;
        myNumTaxa = numTaxa;
        mySiteScoreType = siteScoreType;
        myValues = SuperByteMatrixBuilder.getInstance(myNumTaxa, myNumSites);
        myTaxaList = taxaList;
    }

    public static Byte2DBuilder getInstance(int numTaxa, int numSites, SiteScore.SITE_SCORE_TYPE siteScoreType, TaxaList taxaList) {
        return new Byte2DBuilder(numTaxa, numSites, siteScoreType, taxaList);
    }

    public static Byte2DBuilder getInstance(IHDF5Writer writer, int numSites, SiteScore.SITE_SCORE_TYPE siteScoreType, TaxaList taxaList) {
        return new Byte2DBuilder(writer, numSites, siteScoreType, taxaList);
    }

    public static FilterByte2D getFilteredInstance(Byte2D base, FilterGenotypeTable filterGenotypeTable) {
        return new FilterByte2D(base, filterGenotypeTable);
    }

    public static Byte2D getInstance(IHDF5Reader reader, SiteScore.SITE_SCORE_TYPE siteScoreType) {
        return new HDF5Byte2D(reader, siteScoreType);
    }

    /**
     * Add taxon and set values for all sites for that taxon.
     *
     * @param taxon taxon
     * @param values values
     *
     * @return builder
     */
    public Byte2DBuilder addTaxon(int taxon, byte[] values) {

        if (values.length != myNumSites) {
            throw new IllegalStateException("Byte2DBuilder: addTaxon: Number of sites: " + values.length + " should be: " + myNumSites);
        }

        if (myIsHDF5) {
            synchronized (myHDF5Writer) {
                HDF5Utils.writeHDF5GenotypeSiteScores(myHDF5Writer, myTaxaList.get(taxon).getName(), mySiteScoreType.toString(), values);
            }
            myNumTaxa++;
        } else {
            for (int s = 0; s < myNumSites; s++) {
                myValues.set(taxon, s, values[s]);
            }
        }

        return this;
    }

    public Byte2D build() {
        if (myIsHDF5) {
            IHDF5Reader reader = myHDF5Writer;
            myHDF5Writer = null;
            return new HDF5Byte2D(reader, mySiteScoreType);
        } else {
            SuperByteMatrix temp = myValues;
            myValues = null;
            return new MemoryByte2D(mySiteScoreType, temp);
        }
    }
}
