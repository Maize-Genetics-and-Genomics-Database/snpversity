/*
 *  FilterGenotypeCallTable
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;

/**
 * Filtering class for GenotypeCallTable.  This class is generally never used directly, but rather used through
 * FilterGenotypeTable.
 * @see net.maizegenetics.dna.snp.FilterGenotypeTable
 * @author Terry Casstevens
 */
class FilterGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeCallTable myBaseGenotype;
    private final boolean myIsTaxaFilter;
    private final boolean myIsSiteFilter;
    private final boolean myIsSiteFilterByRange;
    private final int[] myTaxaRedirect;
    private final int[] mySiteRedirect;
    private final int myRangeStart;
    private final int myRangeEnd;

    FilterGenotypeCallTable(GenotypeCallTable genotype, int numTaxa, int[] taxaRedirect, int numSites, int rangeStart, int rangeEnd) {

        super(numTaxa, numSites, genotype.isPhased(), null, genotype.maxNumAlleles());
        myBaseGenotype = genotype;

        myTaxaRedirect = taxaRedirect;
        if (taxaRedirect == null) {
            myIsTaxaFilter = false;
        } else {
            myIsTaxaFilter = true;
        }

        myIsSiteFilter = false;
        mySiteRedirect = null;

        myRangeStart = rangeStart;
        myRangeEnd = rangeEnd;
        if ((rangeStart == -1) || (rangeEnd == -1)) {
            myIsSiteFilterByRange = false;
        } else {
            myIsSiteFilterByRange = true;
        }

    }

    FilterGenotypeCallTable(GenotypeCallTable genotype, int numTaxa, int[] taxaRedirect, int numSites, int[] siteRedirect) {

        super(numTaxa, numSites, genotype.isPhased(), null, genotype.maxNumAlleles());
        myBaseGenotype = genotype;

        myTaxaRedirect = taxaRedirect;
        if (taxaRedirect == null) {
            myIsTaxaFilter = false;
        } else {
            myIsTaxaFilter = true;
        }

        myRangeStart = -1;
        myRangeEnd = -1;
        myIsSiteFilterByRange = false;

        mySiteRedirect = siteRedirect;
        if (siteRedirect == null) {
            myIsSiteFilter = false;
        } else {
            myIsSiteFilter = true;
        }

    }

    @Override
    public byte genotype(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return GenotypeTable.UNKNOWN_ALLELE;
        } else {
            return myBaseGenotype.genotype(taxaIndex, translateSite(site));
        }
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return GenotypeTable.UNKNOWN_ALLELE_STR;
        } else {
            return myBaseGenotype.genotypeAsString(taxaIndex, translateSite(site));
        }
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        if (!myIsTaxaFilter) {
            return myBaseGenotype.genotypeForAllTaxa(translateSite(site));
        } else {
            return super.genotypeForAllTaxa(site);
        }
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String[][] alleleDefinitions() {
        String[][] encodings = myBaseGenotype.alleleDefinitions();
        if (encodings.length == 1) {
            return encodings;
        } else if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = numberOfSites();
            String[][] result = new String[numSites][];
            for (int i = 0; i < numSites; i++) {
                result[i] = alleleDefinitions(i);
            }
            return result;
        } else {
            return encodings;
        }
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return myBaseGenotype.alleleDefinitions(translateSite(site));
    }

    @Override
    public int maxNumAlleles() {
        return myBaseGenotype.maxNumAlleles();
    }

    public int translateSite(int site) {

        if (myIsSiteFilterByRange) {
            return site + myRangeStart;
        } else if (myIsSiteFilter) {
            return mySiteRedirect[site];
        } else {
            return site;
        }

    }

    public int translateTaxon(int taxon) {

        if (myIsTaxaFilter) {
            return myTaxaRedirect[taxon];
        } else {
            return taxon;
        }

    }
}
