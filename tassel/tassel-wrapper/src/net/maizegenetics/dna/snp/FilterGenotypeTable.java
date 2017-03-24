/*
 * FilterGenotypeTable
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.bit.DynamicBitStorage;
import net.maizegenetics.dna.snp.depth.AlleleDepth;
import net.maizegenetics.dna.snp.depth.FilterAlleleDepth;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.AlleleProbabilityBuilder;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder;
import net.maizegenetics.dna.snp.score.Dosage;
import net.maizegenetics.dna.snp.score.DosageBuilder;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.GeneralAnnotationStorage;

import org.apache.log4j.Logger;

import java.util.*;
import java.util.stream.Stream;

/**
 * Taxa and site filtering of GenotypeTables. The class essentially creates
 * views of the baseGenotypeTable through arrays for indirection.
 *
 * @author Terry Casstevens
 */
public class FilterGenotypeTable implements GenotypeTable {

    private static final long serialVersionUID = -5197800047652332969L;
    private static final Logger myLogger = Logger.getLogger(FilterGenotypeTable.class);
    private final boolean myIsTaxaFilter;
    private final boolean myIsSiteFilter;
    private final boolean myIsSiteFilterByRange;
    private final GenotypeTable myBaseAlignment;
    private final TaxaList myTaxaList;
    private final int[] myTaxaRedirect;
    private final int[] mySiteRedirect;
    private final int myRangeStart;
    private final int myRangeEnd;
    private Chromosome[] myChromosomes;
    private int[] myChromosomeOffsets;
    private PositionList myPositionList;
    private AlleleDepth myAlleleDepth = null;
    private AlleleProbability myAlleleProbability;
    private ReferenceProbability myReferenceProbabily;
    private Dosage myDosage;
    private final GenotypeCallTable myGenotype;
    private final Map<WHICH_ALLELE, BitStorage> myBitStorage = new HashMap<>();

    private FilterGenotypeTable(GenotypeTable a, TaxaList subList, int[] taxaRedirect, FilterGenotypeTable original) {

        myTaxaList = subList;

        if (myTaxaList.numberOfTaxa() != taxaRedirect.length) {
            throw new IllegalArgumentException("FilterGenotypeTable: init: subList should be same size as taxaRedirect.");
        }

        myIsTaxaFilter = true;
        myBaseAlignment = a;
        myTaxaRedirect = taxaRedirect;

        if (original == null) {
            myIsSiteFilter = false;
            myIsSiteFilterByRange = false;
            mySiteRedirect = null;
            myRangeStart = -1;
            myRangeEnd = -1;
            myChromosomes = myBaseAlignment.chromosomes();
            myChromosomeOffsets = myBaseAlignment.chromosomesOffsets();
        } else {
            myIsSiteFilter = original.isSiteFilter();
            myIsSiteFilterByRange = original.isSiteFilterByRange();
            mySiteRedirect = original.getSiteRedirect();
            myRangeStart = original.getRangeStart();
            myRangeEnd = original.getRangeEnd();
            myChromosomes = original.chromosomes();
            myChromosomeOffsets = original.chromosomesOffsets();
        }

        GenotypeCallTable myCallTable = myBaseAlignment.genotypeMatrix();
        if (myIsSiteFilter && myCallTable != null) {
            myGenotype = GenotypeCallTableBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), mySiteRedirect);
        } else if (myCallTable != null) {
            myGenotype = GenotypeCallTableBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), myRangeStart, myRangeEnd);
        } else {
            myGenotype = null;
        }

    }

    /**
     * This returns FilterGenotypeTable with only specified subTaxaList.
     * Defaults to retain unknown taxa.
     *
     * @param a alignment
     * @param subTaxaList subset id group
     *
     * @return filter alignment
     */
    public static GenotypeTable getInstance(GenotypeTable a, TaxaList subTaxaList) {
        return getInstance(a, subTaxaList, true);
    }

    /**
     * This returns FilterGenotypeTable with only specified subTaxaList. If
     * retainUnknownTaxa is true then Alignment will return unknown values for
     * missing taxa.
     *
     * @param a alignment
     * @param subTaxaList subset id group
     * @param retainUnknownTaxa whether to retain unknown taxa
     *
     * @return filter alignment
     */
    public static GenotypeTable getInstance(GenotypeTable a, TaxaList subTaxaList, boolean retainUnknownTaxa) {

        GenotypeTable baseAlignment = a;
        FilterGenotypeTable original = null;
        if (baseAlignment instanceof FilterGenotypeTable) {
            original = (FilterGenotypeTable) a;
            baseAlignment = ((FilterGenotypeTable) a).getBaseAlignment();
        }

        List<Integer> taxaRedirectList = new ArrayList<Integer>();
        List<Taxon> idList = new ArrayList<Taxon>();
        boolean noNeedToFilter = true;
        if (subTaxaList.numberOfTaxa() != a.numberOfTaxa()) {
            noNeedToFilter = false;
        }

        for (int i = 0, n = subTaxaList.numberOfTaxa(); i < n; i++) {
            int ion = a.taxa().indexOf(subTaxaList.get(i));

            if (ion != i) {
                noNeedToFilter = false;
            }

            if (ion == -1) {
                if (retainUnknownTaxa) {
                    taxaRedirectList.add(-1);
                    idList.add(subTaxaList.get(i));
                }
            } else {
                if (a instanceof FilterGenotypeTable) {
                    taxaRedirectList.add(((FilterGenotypeTable) a).translateTaxon(ion));
                } else {
                    taxaRedirectList.add(ion);
                }
                idList.add(a.taxa().get(ion));
            }
        }

        if (noNeedToFilter) {
            return a;
        }

        int[] taxaRedirect = new int[taxaRedirectList.size()];
        for (int j = 0, n = taxaRedirectList.size(); j < n; j++) {
            taxaRedirect[j] = (int) taxaRedirectList.get(j);
        }

        TaxaList resultTaxaList = new TaxaListBuilder().addAll(idList).build();

        return new FilterGenotypeTable(baseAlignment, resultTaxaList, taxaRedirect, original);

    }

    /**
     * Removes specified IDs.
     *
     * @param a alignment to filter
     * @param subTaxaList specified IDs
     *
     * @return Filtered Alignment
     */
    public static GenotypeTable getInstanceRemoveIDs(GenotypeTable a, TaxaList subTaxaList) {

        TaxaListBuilder result = new TaxaListBuilder();
        TaxaList current = a.taxa();
        for (int i = 0, n = current.numberOfTaxa(); i < n; i++) {
            if (subTaxaList.indexOf(current.get(i)) < 0) {
                result.add(current.get(i));
            }
        }
        return FilterGenotypeTable.getInstance(a, result.build());

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param startSite start site (included)
     * @param endSite end site (included)
     */
    private FilterGenotypeTable(GenotypeTable a, int startSite, int endSite, FilterGenotypeTable original) {

        myTaxaList = original == null ? a.taxa() : original.taxa();

        if (startSite > endSite) {
            throw new IllegalArgumentException("FilterGenotypeTable: init: start site: " + startSite + " is larger than end site: " + endSite);
        }

        if ((startSite < 0) || (startSite > a.numberOfSites() - 1)) {
            throw new IllegalArgumentException("FilterGenotypeTable: init: start site: " + startSite + " is out of range.");
        }

        if ((endSite < 0) || (endSite > a.numberOfSites() - 1)) {
            throw new IllegalArgumentException("FilterGenotypeTable: init: end site: " + endSite + " is out of range.");
        }

        myBaseAlignment = a;
        myIsSiteFilterByRange = true;
        myIsSiteFilter = false;
        myRangeStart = startSite;
        myRangeEnd = endSite;
        mySiteRedirect = null;
        getLociFromBase();

        if (original == null) {
            myIsTaxaFilter = false;
            myTaxaRedirect = null;
        } else {
            myIsTaxaFilter = original.isTaxaFilter();
            myTaxaRedirect = original.getTaxaRedirect();
        }

        if (myIsSiteFilter) {
            myGenotype = GenotypeCallTableBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), mySiteRedirect);
        } else {
            myGenotype = GenotypeCallTableBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), myRangeStart, myRangeEnd);
        }

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param subSites site to include
     */
    private FilterGenotypeTable(GenotypeTable a, int[] subSites, FilterGenotypeTable original) {

        myTaxaList = original == null ? a.taxa() : original.taxa();

        myBaseAlignment = a;
        myIsSiteFilter = true;
        myIsSiteFilterByRange = false;
        if ((subSites == null) || (subSites.length == 0)) {
            mySiteRedirect = new int[0];
        } else {
            mySiteRedirect = new int[subSites.length];
            Arrays.sort(subSites);
            System.arraycopy(subSites, 0, mySiteRedirect, 0, subSites.length);
        }
        myRangeStart = -1;
        myRangeEnd = -1;
        getLociFromBase();

        if (original == null) {
            myIsTaxaFilter = false;
            myTaxaRedirect = null;
        } else {
            myIsTaxaFilter = original.isTaxaFilter();
            myTaxaRedirect = original.getTaxaRedirect();
        }

        if (myIsSiteFilter) {
            myGenotype = GenotypeCallTableBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), mySiteRedirect);
        } else {
            myGenotype = GenotypeCallTableBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), myRangeStart, myRangeEnd);
        }

    }

    public static FilterGenotypeTable getInstance(GenotypeTable a, int[] subSites) {

        if (a instanceof FilterGenotypeTable) {
            FilterGenotypeTable original = (FilterGenotypeTable) a;
            GenotypeTable baseAlignment = ((FilterGenotypeTable) a).getBaseAlignment();
            if (original.isSiteFilter()) {
                int[] newSubSites = new int[subSites.length];
                for (int i = 0; i < subSites.length; i++) {
                    newSubSites[i] = original.translateSite(subSites[i]);
                }
                return new FilterGenotypeTable(baseAlignment, newSubSites, original);
            } else if (original.isSiteFilterByRange()) {
                int[] newSubSites = new int[subSites.length];
                for (int i = 0; i < subSites.length; i++) {
                    newSubSites[i] = original.translateSite(subSites[i]);
                }
                return new FilterGenotypeTable(baseAlignment, newSubSites, original);
            } else if (original.isTaxaFilter()) {
                return new FilterGenotypeTable(baseAlignment, subSites, original);
            } else {
                throw new IllegalStateException("FilterGenotypeTable: getInstance: original not in known state.");
            }
        } else {
            return new FilterGenotypeTable(a, subSites, null);
        }

    }

    public static FilterGenotypeTable getInstance(GenotypeTable a, List<String> siteNamesToKeep) {
        return getInstance(a, siteNamesToKeep.toArray(new String[siteNamesToKeep.size()]));
    }

    public static FilterGenotypeTable getInstance(GenotypeTable a, String[] siteNamesToKeep) {

        Arrays.sort(siteNamesToKeep);
        int[] temp = new int[siteNamesToKeep.length];
        int count = 0;
        for (int i = 0, n = a.numberOfSites(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToKeep, a.siteName(i)) >= 0) {
                temp[count++] = i;
                if (count == siteNamesToKeep.length) {
                    break;
                }
            }
        }

        int[] result = null;
        if (count == siteNamesToKeep.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static FilterGenotypeTable getInstanceRemoveSiteNames(GenotypeTable a, List<String> siteNamesToRemove) {
        return getInstance(a, siteNamesToRemove.toArray(new String[siteNamesToRemove.size()]));
    }

    public static FilterGenotypeTable getInstanceRemoveSiteNames(GenotypeTable a, String[] siteNamesToRemove) {

        Arrays.sort(siteNamesToRemove);
        int[] temp = new int[a.numberOfSites()];
        int count = 0;
        for (int i = 0, n = a.numberOfSites(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToRemove, a.siteName(i)) < 0) {
                temp[count++] = i;
            }
        }

        int[] result = null;
        if (count == temp.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static FilterGenotypeTable getInstance(GenotypeTable a, PositionList subPositionList) {

        int[] temp = new int[subPositionList.size()];
        int count = 0;
        PositionList positionList = a.positions();
        for (Position position : subPositionList) {
            int index = positionList.indexOf(position);
            if (index >= 0) {
                temp[count++] = index;
            }
        }

        int[] result = null;
        if (count == subPositionList.size()) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static FilterGenotypeTable getInstance(GenotypeTable a, String chromosome, int startPhysicalPos, int endPhysicalPos) {
        return getInstance(a, a.chromosome(chromosome), startPhysicalPos, endPhysicalPos);
    }

    public static FilterGenotypeTable getInstance(GenotypeTable a, Chromosome chromosome, int startPhysicalPos, int endPhysicalPos) {

        int startSite = a.siteOfPhysicalPosition(startPhysicalPos, chromosome);
        if (startSite < 0) {
            startSite = -(startSite + 1);
        }

        int endSite = a.siteOfPhysicalPosition(endPhysicalPos, chromosome);
        if (endSite < 0) {
            endSite = -(endSite + 2);
        }

        if (startSite > endSite) {
            myLogger.warn("getInstance: start site: " + startSite + " from physical pos: " + startPhysicalPos + " is larger than end site: " + endSite + " from physical pos: " + endPhysicalPos);
            return null;
        }

        return getInstance(a, startSite, endSite);

    }

    public static FilterGenotypeTable getInstance(GenotypeTable a, Chromosome chromosome) {
        int[] endStart = a.firstLastSiteOfChromosome(chromosome);
        return getInstance(a, endStart[0], endStart[1]);
    }

    /**
     * Factory method that returns a FilterGenotypeTable viewing sites between
     * start site and end site inclusive.
     *
     * @param a alignment
     * @param startSite start site
     * @param endSite end site
     *
     * @return Filter Alignment
     */
    public static FilterGenotypeTable getInstance(GenotypeTable a, int startSite, int endSite) {

        if (a instanceof FilterGenotypeTable) {
            FilterGenotypeTable original = (FilterGenotypeTable) a;
            GenotypeTable baseAlignment = ((FilterGenotypeTable) a).getBaseAlignment();
            if (original.isSiteFilter()) {
                int[] subSites = new int[endSite - startSite + 1];
                int[] originalSites = original.getSiteRedirect();
                for (int i = startSite; i <= endSite; i++) {
                    subSites[i - startSite] = originalSites[i];
                }
                return new FilterGenotypeTable(baseAlignment, subSites, original);
            } else if (original.isSiteFilterByRange()) {
                return new FilterGenotypeTable(baseAlignment, original.translateSite(startSite), original.translateSite(endSite), original);
            } else if (original.isTaxaFilter()) {
                return new FilterGenotypeTable(baseAlignment, startSite, endSite, original);
            } else {
                throw new IllegalStateException("FilterGenotypeTable: getInstance: original not in known state.");
            }
        } else {
            return new FilterGenotypeTable(a, startSite, endSite, null);
        }

    }

    @Override
    public byte genotype(int taxon, int site) {
        return myGenotype.genotype(taxon, site);
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeRange(taxon, startSite, endSite);
    }

    @Override
    public byte genotype(int taxon, Chromosome chromosome, int physicalPosition) {
        return myGenotype.genotype(taxon, myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome));
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

    /**
     * Returns site of this FilterGenotypeTable based on given site from
     * embedded Alignment.
     *
     * @param site site
     * @return site in this alignment
     */
    public int reverseTranslateSite(int site) {

        if (myIsSiteFilterByRange) {
            return site - myRangeStart;
        } else if (myIsSiteFilter) {
            return Arrays.binarySearch(mySiteRedirect, site);
        } else {
            return site;
        }

    }

    /**
     * Returns sites from original alignment that are viewable (not filtered) by
     * this filter alignment.
     *
     * @return list of sites
     */
    public int[] getBaseSitesShown() {
        int numSites = numberOfSites();
        int[] result = new int[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = translateSite(i);
        }
        return result;
    }

    public int translateTaxon(int taxon) {

        if (myIsTaxaFilter) {
            return myTaxaRedirect[taxon];
        } else {
            return taxon;
        }

    }

    private void getLociFromBase() {

        if ((!myIsSiteFilter) && (!myIsSiteFilterByRange)) {
            myChromosomes = myBaseAlignment.chromosomes();
            myChromosomeOffsets = myBaseAlignment.chromosomesOffsets();
            return;
        }

        int numSites = numberOfSites();
        List<Chromosome> chromosomes = new ArrayList<Chromosome>();
        List<Integer> offsets = new ArrayList<Integer>();
        for (int i = 0; i < numSites; i++) {
            Chromosome current = chromosome(i);
            if (!chromosomes.contains(current)) {
                chromosomes.add(current);
                offsets.add(i);
            }
        }

        myChromosomes = new Chromosome[chromosomes.size()];
        chromosomes.toArray(myChromosomes);

        myChromosomeOffsets = new int[offsets.size()];
        for (int i = 0; i < offsets.size(); i++) {
            myChromosomeOffsets[i] = (Integer) offsets.get(i);
        }

    }

    @Override
    public int indelSize(int site) {
        return myBaseAlignment.indelSize(translateSite(site));
    }

    @Override
    public Chromosome chromosome(int site) {
        return myBaseAlignment.chromosome(translateSite(site));
    }

    @Override
    public int chromosomalPosition(int site) {
        return myBaseAlignment.chromosomalPosition(translateSite(site));
    }

    @Override
    public String chromosomeName(int site) {
        return myBaseAlignment.chromosomeName(translateSite(site));
    }

    @Override
    public Chromosome chromosome(String name) {
        return myBaseAlignment.chromosome(name);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myChromosomes;
    }

    @Override
    public int numChromosomes() {
        return myChromosomes.length;
    }

    @Override
    public int[] chromosomesOffsets() {
        return myChromosomeOffsets;
    }

    @Override
    public int[] firstLastSiteOfChromosome(Chromosome chromosome) {
        for (int i = 0; i < numChromosomes(); i++) {
            if (chromosome.equals(myChromosomes[i])) {
                int end = 0;
                if (i == numChromosomes() - 1) {
                    end = numberOfSites() - 1;
                } else {
                    end = myChromosomeOffsets[i + 1] - 1;
                }
                return new int[]{myChromosomeOffsets[i], end};
            }
        }
        throw new IllegalArgumentException("FilterGenotypeTable: getStartAndEndOfLocus: this locus not defined: " + chromosome.getName());
    }

    @Override
    public byte referenceAllele(int site) {
        return myBaseAlignment.referenceAllele(translateSite(site));
    }

    @Override
    public int numberOfSites() {

        if (myIsSiteFilterByRange) {
            return myRangeEnd - myRangeStart + 1;
        } else if (myIsSiteFilter) {
            return mySiteRedirect.length;
        } else {
            return myBaseAlignment.numberOfSites();
        }

    }

    @Override
    public String siteName(int site) {
        return myBaseAlignment.siteName(translateSite(site));
    }

    @Override
    public boolean hasReference() {
        return myBaseAlignment.hasReference();
    }

    @Override
    public boolean isIndel(int site) {
        return myBaseAlignment.isIndel(translateSite(site));
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        int temp = myBaseAlignment.siteOfPhysicalPosition(physicalPosition, chromosome);
        if (temp < 0) {
            temp = -(temp + 1);
            return -(reverseTranslateSite(temp) + 1);
        }
        return reverseTranslateSite(temp);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        int temp = myBaseAlignment.siteOfPhysicalPosition(physicalPosition, chromosome, snpName);
        if (temp < 0) {
            temp = -(temp + 1);
            return -(reverseTranslateSite(temp) + 1);
        }
        return reverseTranslateSite(temp);
    }

    public GenotypeTable getBaseAlignment() {
        return myBaseAlignment;
    }

    public boolean isTaxaFilter() {
        return myIsTaxaFilter;
    }

    public boolean isSiteFilter() {
        return myIsSiteFilter;
    }

    public boolean isSiteFilterByRange() {
        return myIsSiteFilterByRange;
    }

    protected int[] getTaxaRedirect() {
        return myTaxaRedirect;
    }

    protected int[] getSiteRedirect() {
        return mySiteRedirect;
    }

    protected int getRangeStart() {
        return myRangeStart;
    }

    protected int getRangeEnd() {
        return myRangeEnd;
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return myGenotype.genotypeArray(taxon, site);
    }

    @Override
    public byte[] genotypeAllTaxa(int site) {
        return myGenotype.genotypeForAllTaxa(site);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        return myGenotype.genotypeAllSites(taxon);
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, WHICH_ALLELE allele) {
        return bitStorage(allele).allelePresenceForAllSites(taxon);
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, WHICH_ALLELE allele, int startBlock, int endBlock) {
        return bitStorage(allele).allelePresenceForSitesBlock(taxon, startBlock, endBlock);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myGenotype.genotypeAsString(taxon, site);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        return myGenotype.genotypeAsStringArray(taxon, site);
    }

    @Override
    public byte[] referenceAlleleForAllSites() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            byte[] result = new byte[numberOfSites()];
            for (int i = 0, n = numberOfSites(); i < n; i++) {
                result[i] = referenceAllele(i);
            }
            return result;
        } else {
            return myBaseAlignment.referenceAlleleForAllSites();
        }
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        return myGenotype.isHeterozygous(taxon, site);
    }

    @Override
    public int[] physicalPositions() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = numberOfSites();
            int[] result = new int[numSites];
            for (int i = 0; i < numSites; i++) {
                result[i] = chromosomalPosition(i);
            }
            return result;
        } else {
            return myBaseAlignment.physicalPositions();
        }
    }

    @Override
    public TaxaList taxa() {
        return myTaxaList;
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaList.numberOfTaxa();
    }

    @Override
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        return myBaseAlignment.siteScoreTypes();
    }

    @Override
    public String genomeVersion() {
        return myBaseAlignment.genomeVersion();
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myBaseAlignment.isPositiveStrand(translateSite(site));
    }

    @Override
    public boolean isPhased() {
        return myGenotype.isPhased();
    }

    @Override
    public boolean retainsRareAlleles() {
        return myGenotype.retainsRareAlleles();
    }

    @Override
    public String[][] alleleDefinitions() {
        return myGenotype.alleleDefinitions();
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return myGenotype.alleleDefinitions(site);
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        return myGenotype.genotypeAsString(site, value);
    }

    @Override
    public int maxNumAlleles() {
        return myGenotype.maxNumAlleles();
    }

    @Override
    public byte[] alleles(int site) {
        return myGenotype.alleles(site);
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return myGenotype.allelesSortedByFrequency(site);
    }

    @Override
    public double majorAlleleFrequency(int site) {
        return myGenotype.majorAlleleFrequency(site);
    }

    @Override
    public int heterozygousCount(int site) {
        return myGenotype.heterozygousCount(site);
    }

    @Override
    public boolean isPolymorphic(int site) {
        return myGenotype.isPolymorphic(site);
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        return myGenotype.totalGametesNonMissingForSite(site);
    }

    @Override
    public int totalNonMissingForSite(int site) {
        return myGenotype.totalNonMissingForSite(site);
    }

    @Override
    public int minorAlleleCount(int site) {
        return myGenotype.minorAlleleCount(site);
    }

    @Override
    public double minorAlleleFrequency(int site) {
        return myGenotype.minorAlleleFrequency(site);
    }

    @Override
    public int majorAlleleCount(int site) {
        return myGenotype.majorAlleleCount(site);
    }

    @Override
    public byte majorAllele(int site) {
        return myGenotype.majorAllele(site);
    }

    @Override
    public byte minorAllele(int site) {
        return myGenotype.minorAllele(site);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        return myGenotype.genosSortedByFrequency(site);
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        return myGenotype.totalGametesNonMissingForTaxon(taxon);
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        return myGenotype.totalNonMissingForTaxon(taxon);
    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        return myGenotype.heterozygousCountForTaxon(taxon);
    }

    @Override
    public boolean hasGenotype() {
        return myBaseAlignment.hasGenotype();
    }

    @Override
    public boolean hasDepth() {
        return myBaseAlignment.hasDepth();
    }

    @Override
    public boolean hasAlleleProbabilities() {
        return myBaseAlignment.hasAlleleProbabilities();
    }

    @Override
    public boolean hasReferenceProbablity() {
        return myBaseAlignment.hasReferenceProbablity();
    }

    @Override
    public boolean hasDosage() {
        return myBaseAlignment.hasDosage();
    }

    @Override
    public AlleleDepth depth() {
        AlleleDepth depth = myBaseAlignment.depth();
        if (depth == null) {
            return null;
        } else if (myAlleleDepth == null) {
            myAlleleDepth = new FilterAlleleDepth(depth, this);
        }
        return myAlleleDepth;
    }

    @Override
    public int[] depthForAlleles(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return null;
        } else {
            return myBaseAlignment.depthForAlleles(taxaIndex, translateSite(site));
        }
    }

    @Override
    public byte[] allelesBySortType(ALLELE_SORT_TYPE scope, int site) {
        if (scope == ALLELE_SORT_TYPE.Frequency) {
            return alleles(site);
        } else {
            return myBaseAlignment.allelesBySortType(scope, translateSite(site));
        }
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, WHICH_ALLELE allele) {
        return bitStorage(allele).allelePresenceForAllTaxa(site);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, WHICH_ALLELE allele) {
        return bitStorage(allele).haplotypeAllelePresenceForAllSites(taxon, firstParent);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, WHICH_ALLELE allele) {
        return bitStorage(allele).haplotypeAllelePresenceForAllTaxa(site, firstParent);
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, WHICH_ALLELE allele, int startBlock, int endBlock) {
        return bitStorage(allele).haplotypeAllelePresenceForSitesBlock(taxon, firstParent, startBlock, endBlock);
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeAsStringRange(taxon, startSite, endSite);
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        return myGenotype.genotypeAsStringRow(taxon);
    }

    @Override
    public byte[] referenceAlleles(int startSite, int endSite) {
        if (!hasReference()) {
            return null;
        }

        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i] = referenceAllele(i);
        }
        return result;
    }

    @Override
    public int chromosomeSiteCount(Chromosome chromosome) {
        int[] startEnd = firstLastSiteOfChromosome(chromosome);
        return startEnd[1] - startEnd[0] + 1;
    }

    @Override
    public boolean isAllPolymorphic() {
        return myGenotype.isAllPolymorphic();
    }

    @Override
    public String majorAlleleAsString(int site) {
        return myGenotype.majorAlleleAsString(site);
    }

    @Override
    public String minorAlleleAsString(int site) {
        return myGenotype.minorAlleleAsString(site);
    }

    @Override
    public byte[] minorAlleles(int site) {
        return myGenotype.minorAlleles(site);
    }

    @Override
    public String taxaName(int index) {
        return myTaxaList.taxaName(index);
    }

    @Override
    public GenotypeTable[] compositeAlignments() {
        return new GenotypeTable[]{this};
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myGenotype.diploidAsString(site, value);
    }

    @Override
    public Object[][] genoCounts() {
        return myGenotype.genoCounts();
    }

    @Override
    public Object[][] majorMinorCounts() {
        return myGenotype.majorMinorCounts();
    }

    @Override
    public BitStorage bitStorage(WHICH_ALLELE allele) {
        BitStorage result = myBitStorage.get(allele);
        if (result != null) {
            return result;
        }

        switch (allele) {
            case Major:
                result = new DynamicBitStorage(myGenotype, allele, myGenotype.majorAlleleForAllSites());
                break;
            case Minor:
                result = new DynamicBitStorage(myGenotype, allele, myGenotype.minorAlleleForAllSites());
                break;
            case Minor2:
                result = new DynamicBitStorage(myGenotype, allele, myGenotype.thirdAlleleForAllSites());
                break;
            case Unknown:
                result = DynamicBitStorage.getUnknownInstance(myGenotype);
                break;
            default:
                myLogger.warn("bitStorage: Unsupported allele: " + allele);
                return null;
        }

        myBitStorage.put(allele, result);
        return result;
    }

    @Override
    public PositionList positions() {
        if (myPositionList == null) {
            PositionListBuilder pLB = new PositionListBuilder();
            PositionList basePL = getBaseAlignment().positions();
            for (int i = 0; i < numberOfSites(); i++) {
                pLB.add(basePL.get(translateSite(i)));
            }
            myPositionList = pLB.build();
        }
        return myPositionList;
    }

    @Override
    public GenotypeCallTable genotypeMatrix() {
        return myGenotype;
    }

    @Override
    public AlleleProbability alleleProbability() {
        AlleleProbability probability = myBaseAlignment.alleleProbability();
        if (probability == null) {
            return null;
        } else if (myAlleleProbability == null) {
            myAlleleProbability = AlleleProbabilityBuilder.getFilteredInstance(probability, this);
        }
        return myAlleleProbability;
    }

    @Override
    public float alleleProbability(int taxon, int site, SITE_SCORE_TYPE type) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Float.NaN;
        } else {
            return myBaseAlignment.alleleProbability(taxaIndex, translateSite(site), type);
        }
    }

    @Override
    public ReferenceProbability referenceProbability() {
        ReferenceProbability probability = myBaseAlignment.referenceProbability();
        if (probability == null) {
            return null;
        } else if (myReferenceProbabily == null) {
            myReferenceProbabily = ReferenceProbabilityBuilder.getFilteredInstance(probability, this);
        }
        return myReferenceProbabily;
    }

    @Override
    public float referenceProbability(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Float.NaN;
        } else {
            return myBaseAlignment.referenceProbability(taxaIndex, translateSite(site));
        }
    }

    @Override
    public Dosage dosage() {
        Dosage dosage = myBaseAlignment.dosage();
        if (dosage == null) {
            return null;
        } else if (myDosage == null) {
            myDosage = DosageBuilder.getFilteredInstance(dosage, this);
        }
        return myDosage;
    }

    @Override
    public byte dosage(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return 0;
        } else {
            return myBaseAlignment.dosage(taxaIndex, translateSite(site));
        }
    }

    @Override
    public GeneralAnnotationStorage annotations() {
        return myBaseAlignment.annotations();
    }

    @Override
    public Stream<Byte> streamGenotype() {
        return myGenotype.stream();
    }

    @Override
    public Stream<Byte> streamGenotype(int taxon) {
        return myGenotype.stream(taxon);
    }
    
}
