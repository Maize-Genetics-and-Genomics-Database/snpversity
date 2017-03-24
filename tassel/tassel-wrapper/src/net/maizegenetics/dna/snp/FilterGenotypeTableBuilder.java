/*
 *  FilterGenotypeTableBuilder
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;

import java.util.List;

/**
 *
 * @author Terry Casstevens
 */
public class FilterGenotypeTableBuilder {

    private final GenotypeTable myBaseGenotypeTable;

    // Site Related
    private double myMinFreqForSite = 0.0;
    private double myMaxFreqForSite = 1.0;
    private int myMinCountForSite = 0;
    private int[] mySitesToKeep = null;
    private List<String> mySiteNamesToKeep = null;
    private List<String> mySiteNamesToRemove = null;

    // Taxa Related
    private TaxaList myTaxaToKeep = null;
    private TaxaList myTaxaToRemove = null;
    private double myMinHeterozygousForTaxon = 0.0;
    private double myMaxHeterozygousForTaxon = 1.0;
    private int myMinNotMissingForTaxon = 0;

    private FilterGenotypeTableBuilder(GenotypeTable baseGenotypeTable) {
        myBaseGenotypeTable = baseGenotypeTable;
    }

    /**
     * @deprecated
     * 
     * This class to be replaced by FilterSiteBuilderPlugin and
     * FilterTaxaBuilderPlugin
     * 
     */
    public static FilterGenotypeTableBuilder getInstance(GenotypeTable baseGenotypeTable) {
        return new FilterGenotypeTableBuilder(baseGenotypeTable);
    }

    /**
     * Set minimum and maximum minor allele frequency that site must have to
     * remain after filtering.
     *
     * @param minFreq minimum minor allele frequency (default 0.0)
     * @param maxFreq maximum minor allele frequency (default 1.0)
     *
     * @return this builder
     */
    public FilterGenotypeTableBuilder minorAlleleFreqForSite(double minFreq, double maxFreq) {
        if ((minFreq < 0.0) || (minFreq > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: minorAlleleFreqForSite: Min Value must be between 0.0 and 1.0: " + minFreq);
        }
        if ((maxFreq < 0.0) || (maxFreq > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: minorAlleleFreqForSite: Max Value must be between 0.0 and 1.0: " + maxFreq);
        }
        myMinFreqForSite = minFreq;
        myMaxFreqForSite = maxFreq;
        return this;
    }

    /**
     * Set minimum number of taxa with a known base (i.e. not Unknown (N)).
     *
     * @param count count (default 0)
     *
     * @return this builder
     */
    public FilterGenotypeTableBuilder minCountForSite(int count) {
        myMinCountForSite = count;
        return this;
    }

    public FilterGenotypeTableBuilder sitesToKeep(int[] sitesToKeep) {
        mySitesToKeep = sitesToKeep;
        return this;
    }

    public void siteNamesToKeep(List<String> sitesToKeep) {
        mySiteNamesToKeep = sitesToKeep;
    }

    public void siteNamesToRemove(List<String> sitesToRemove) {
        mySiteNamesToRemove = sitesToRemove;
    }

    public FilterGenotypeTableBuilder taxaToKeep(TaxaList taxaToKeep) {
        myTaxaToKeep = taxaToKeep;
        return this;
    }

    public FilterGenotypeTableBuilder taxaToRemove(TaxaList taxaToRemove) {
        myTaxaToRemove = taxaToRemove;
        return this;
    }

    public FilterGenotypeTableBuilder minNotMissingForTaxon(int minNotMissing) {
        if ((minNotMissing < 0.0) || (minNotMissing > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMinNotMissingForTaxon: Value must be between 0.0 and 1.0: " + minNotMissing);
        }
        myMinNotMissingForTaxon = minNotMissing;
        return this;
    }

    public FilterGenotypeTableBuilder minHeterozygousForTaxon(double minHeterozygous) {
        if ((minHeterozygous < 0.0) || (minHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMinHeterozygousForTaxon: Value must be between 0.0 and 1.0: " + minHeterozygous);
        }
        myMinHeterozygousForTaxon = minHeterozygous;
        return this;
    }

    public FilterGenotypeTableBuilder maxHeterozygousForTaxon(double maxHeterozygous) {
        if ((maxHeterozygous < 0.0) || (maxHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMaxHeterozygousForTaxon: Value must be between 0.0 and 1.0: " + maxHeterozygous);
        }
        myMaxHeterozygousForTaxon = maxHeterozygous;
        return this;
    }

    private GenotypeTable getFilteredAlignment(GenotypeTable genotypeTable) {
        int numSites = genotypeTable.numberOfSites();
        int numTaxa = genotypeTable.numberOfTaxa();
        TaxaList ids = genotypeTable.taxa();

        TaxaListBuilder keepTaxaList = new TaxaListBuilder();
        for (int t = 0; t < numTaxa; t++) {

            if (myMinNotMissingForTaxon != 0.0) {
                int totalNotMissing = genotypeTable.totalNonMissingForTaxon(t);
                double percentNotMissing = (double) totalNotMissing / (double) numSites;
                if (percentNotMissing < myMinNotMissingForTaxon) {
                    continue;
                }
            }

            if ((myMinHeterozygousForTaxon != 0.0) || (myMaxHeterozygousForTaxon != 1.0)) {
                int numHeterozygous = genotypeTable.heterozygousCountForTaxon(t);
                int totalSitesNotMissing = genotypeTable.totalNonMissingForTaxon(t);
                double percentHets = (double) numHeterozygous / (double) totalSitesNotMissing;
                if ((percentHets < myMinHeterozygousForTaxon) || (percentHets > myMaxHeterozygousForTaxon)) {
                    continue;
                }
            }

            keepTaxaList.add(ids.get(t));
        }
        return FilterGenotypeTable.getInstance(genotypeTable, keepTaxaList.build(), false);
    }

    public GenotypeTable build() {

        GenotypeTable result = myBaseGenotypeTable;

        if ((myTaxaToKeep != null) && (myTaxaToRemove != null)) {
            throw new IllegalStateException("FilterGenotypeTableBuilder: build: Taxa list to keep and list to remove defined.");
        } else if (myTaxaToKeep != null) {
            result = FilterGenotypeTable.getInstance(result, myTaxaToKeep, false);
        } else if (myTaxaToRemove != null) {
            result = FilterGenotypeTable.getInstanceRemoveIDs(result, myTaxaToRemove);
        }

        if ((myMinHeterozygousForTaxon != 0.0) || (myMaxHeterozygousForTaxon != 1.0) || (myMinNotMissingForTaxon != 0)) {
            result = getFilteredAlignment(result);
        }

        int count = 0;
        if ((mySitesToKeep != null) && (mySitesToKeep.length != 0)) {
            count++;
        }
        if ((mySiteNamesToKeep != null) && (!mySiteNamesToKeep.isEmpty())) {
            count++;
        }
        if ((mySiteNamesToRemove != null) && (!mySiteNamesToRemove.isEmpty())) {
            count++;
        }

        if (count > 1) {
            throw new IllegalStateException("FilterGenotypeTableBuilder: build: Can only set one of the following: sites to keep, site names to keep, or site names to remove.");
        }

        if (((mySitesToKeep != null) && (mySitesToKeep.length != 0))) {
            result = FilterGenotypeTable.getInstance(result, mySitesToKeep);
        } else if (((mySiteNamesToKeep != null) && (!mySiteNamesToKeep.isEmpty()))) {
            result = FilterGenotypeTable.getInstance(result, mySiteNamesToKeep);
        } else if (((mySiteNamesToRemove != null) && (!mySiteNamesToRemove.isEmpty()))) {
            result = FilterGenotypeTable.getInstanceRemoveSiteNames(result, mySiteNamesToRemove);
        }

        if ((myMinFreqForSite != 0.0) || (myMaxFreqForSite != 1.0) || (myMinCountForSite != 0)) {
            int[] sitesToInclude = GenotypeTableUtils.getIncludedSitesBasedOnFreqIgnoreMissing(result, myMinFreqForSite, myMaxFreqForSite, myMinCountForSite);
            result = FilterGenotypeTable.getInstance(result, sitesToInclude);
        }

        return result;

    }

}
