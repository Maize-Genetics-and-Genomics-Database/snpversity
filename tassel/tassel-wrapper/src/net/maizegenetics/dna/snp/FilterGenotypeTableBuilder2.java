/*
 *  FilterGenotypeTableBuilder
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Terry Casstevens; Jason Wallace
 */
//TODO: Check for null cases of names, numbers, etc, and add in error handling for them
//TODO: Make sure get/set functions are all implemented and logical
//TODO: Wrap in a plugin so can be accessed from CLI and GUI
public class FilterGenotypeTableBuilder2 {

    private final GenotypeTable myBaseGenotypeTable;

    // Site Related
    private double myMinFreqForSite = 0.0;
    private double myMaxFreqForSite = 1.0;
    private int myMinCountForSite = 0;
    private int myMaxCountForSite = Integer.MAX_VALUE;
    private double myMinHeterozygousForSite = 0.0;
    private double myMaxHeterozygousForSite = 1.0;
    private int[] mySitesToKeep = null;
    private List<String> mySiteNamesToKeep = null;
    private List<String> mySiteNamesToRemove = null;

    // Taxa Related
    private TaxaListBuilder myTaxaToKeep = new TaxaListBuilder();
    private TaxaListBuilder myTaxaToRemove = new TaxaListBuilder();
    private double myMinHeterozygousForTaxon = 0.0;
    private double myMaxHeterozygousForTaxon = 1.0;
    private double myMinNotMissingForTaxon = 0.0;
    private double myMaxNotMissingForTaxon = 1.0;

    private FilterGenotypeTableBuilder2(GenotypeTable baseGenotypeTable) {
        myBaseGenotypeTable = baseGenotypeTable;
    }


    public static FilterGenotypeTableBuilder2 getInstance(GenotypeTable baseGenotypeTable) {
        return new FilterGenotypeTableBuilder2(baseGenotypeTable);
    }

    /**
     * Set minimum and maximum minor allele frequency that site must have to
     * remain after filtering.
     *
     * @param minFreq minimum minor allele frequency (default 0.0)
     * @param maxFreq maximum minor allele frequency (default 1.0)
     * @return this builder
     */
    public FilterGenotypeTableBuilder2 minorAlleleFreqForSite(double minFreq, double maxFreq) {
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
     * @return this builder
     */
    public FilterGenotypeTableBuilder2 minCountForSite(int count) {
        myMinCountForSite = count;
        return this;
    }

    public void sitesToKeep(int[] sitesToKeep) {
        mySitesToKeep = sitesToKeep;
    }

    public void siteNamesToKeep(List<String> sitesToKeep) {
        mySiteNamesToKeep = sitesToKeep;
    }

    public void siteNamesToRemove(List<String> sitesToRemove) {
        mySiteNamesToRemove = sitesToRemove;
    }

    public FilterGenotypeTableBuilder2 taxaToKeep(TaxaList taxaToKeep) {
        myTaxaToKeep.addAll(taxaToKeep);
        return this;
    }

    public FilterGenotypeTableBuilder2 taxaToRemove(TaxaList taxaToRemove) {
        myTaxaToRemove.addAll(taxaToRemove);
        return this;
    }

    public FilterGenotypeTableBuilder2 minNotMissingForTaxon(int minNotMissing) {
        if ((minNotMissing < 0.0) || (minNotMissing > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMinNotMissingForTaxon: Value must be between 0.0 and 1.0: " + minNotMissing);
        }
        myMinNotMissingForTaxon = minNotMissing;
        return this;
    }

    public FilterGenotypeTableBuilder2 minHeterozygousForTaxon(double minHeterozygous) {
        if ((minHeterozygous < 0.0) || (minHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMinHeterozygousForTaxon: Value must be between 0.0 and 1.0: " + minHeterozygous);
        }
        myMinHeterozygousForTaxon = minHeterozygous;
        return this;
    }

    public FilterGenotypeTableBuilder2 maxHeterozygousForTaxon(double maxHeterozygous) {
        if ((maxHeterozygous < 0.0) || (maxHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMaxHeterozygousForTaxon: Value must be between 0.0 and 1.0: " + maxHeterozygous);
        }
        myMaxHeterozygousForTaxon = maxHeterozygous;
        return this;
    }

    public FilterGenotypeTableBuilder2 minHeterozygousForSite(double maxHeterozygous) {
        if ((maxHeterozygous < 0.0) || (maxHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMaxHeterozygousForTaxon: Value must be between 0.0 and 1.0: " + maxHeterozygous);
        }
        myMinHeterozygousForSite = maxHeterozygous;
        return this;
    }

    public FilterGenotypeTableBuilder2 maxHeterozygousForSite(double maxHeterozygous) {
        if ((maxHeterozygous < 0.0) || (maxHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterGenotypeTableBuilder: setMaxHeterozygousForTaxon: Value must be between 0.0 and 1.0: " + maxHeterozygous);
        }
        myMaxHeterozygousForSite = maxHeterozygous;
        return this;
    }

    /*private GenotypeTable getFilteredAlignment(GenotypeTable genotypeTable) {
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
    }*/

    public GenotypeTable build() {

        correctNullValues();    //TODO: Add unit test for what happens if null values are passed

        //TODO: UNIT TEST THE HECK OUT OF THIS THING
        TaxaList taxaToKeep = getTaxaListToKeep();
        PositionList sitesToKeep = getPositionListToKeep();

        //TODO: Check that these work as intended, esp. that TaxaList and PositionList have set-like behavior
        GenotypeTable result = myBaseGenotypeTable;
        if(taxaToKeep.size() != result.numberOfTaxa()){
            result = FilterGenotypeTable.getInstance(result, taxaToKeep);
        }
        if(sitesToKeep.size() != result.numberOfSites()){
            result = FilterGenotypeTable.getInstance(result, sitesToKeep);
        }
        return result;
    }

    /**
     * Parse the parameters passed to this plugin and return a TaxaList containing the Taxa to keep
     * @return A TaxaList fo the taxa that passed filtering
     */
    private TaxaList getTaxaListToKeep(){

        //First test if any taxa are specified by name to both keep AND remove; throw error if any found
        TaxaList goodTaxaNames = myTaxaToKeep.build();
        TaxaList badTaxaNames = myTaxaToRemove.build();
        TaxaList taxaConflicts = taxaListsOverlap(goodTaxaNames, badTaxaNames);
        if (taxaConflicts.size() > 0) {
            StringBuilder intersects = new StringBuilder();
            for (Taxon t : taxaConflicts) {
                intersects.append(t.getName() + "\n");
            }
            throw new IllegalStateException("FilterGenotypeTableBuilder: build: Taxa list to keep and to remove contain the same taxa:\n" + intersects.toString());
        }

        //Next build a temporary holder for all the taxa to remove.
        //    (since there's no Builder.remove() function, it's easier to accumulate the ones that fail)
        TaxaListBuilder badTaxaBuilder = new TaxaListBuilder();

        for(int t=0; t<myBaseGenotypeTable.numberOfTaxa(); t++){
            Taxon currentTaxon = myBaseGenotypeTable.taxa().get(t);
            int totalSitesNotMissing = myBaseGenotypeTable.totalNonMissingForTaxon(t);

            //Missingness
            //TODO: Add support for both fractions and integers
            double percentNotMissing = (double) totalSitesNotMissing / (double) myBaseGenotypeTable.numberOfSites();
            if (valueOutsideOfRange(percentNotMissing, myMinNotMissingForTaxon, myMaxNotMissingForTaxon)) {
                badTaxaBuilder.add(currentTaxon);
                continue;
            }

            //Heterozygosity
            int numHeterozygous = myBaseGenotypeTable.heterozygousCountForTaxon(t);
            double percentHets = (double) numHeterozygous / (double) totalSitesNotMissing;
            if ((percentHets < myMinHeterozygousForTaxon) || (percentHets > myMaxHeterozygousForTaxon)) {
                badTaxaBuilder.add(currentTaxon);
                continue;
            }
        }

        //Add in the ones specified by name, then build
        badTaxaBuilder.addAll(badTaxaNames);
        TaxaList badTaxa = badTaxaBuilder.build();

        //Go through and invert the list, returning all (good) taxa that are EITHER specified to keep by name OR are not on the bad list
        TaxaListBuilder goodTaxaBuilder = new TaxaListBuilder();
        for(Taxon t: myBaseGenotypeTable.taxa()){
            if(goodTaxaNames.contains(t) || !badTaxa.contains(t)){
                goodTaxaBuilder.add(t);
            }
        }
        return goodTaxaBuilder.build();
    }


    private void correctNullValues(){
        if(mySitesToKeep == null){
            mySitesToKeep = new int[] {};
        }
        if(mySiteNamesToKeep == null){
            mySiteNamesToKeep = new ArrayList();
        }
        if(mySiteNamesToRemove == null){
            mySiteNamesToRemove = new ArrayList();
        }
    }

    /**
     * Parse the parameters passed to this plugin and return a PositionList of which positions to keep
     * @return A PositionList of sites that passed filtering
     */
    private PositionList getPositionListToKeep(){

        //First test if any sites are specified by name to both keep AND remove; throw error if any found
        //Consolidate the manually specified ones to keep
        PositionList goodPositionsFromNames = getPositionsFromSiteNames(mySiteNamesToKeep, myBaseGenotypeTable);
        PositionList goodPositionsFromIndex = getPositionsFromSiteIndex(mySitesToKeep, myBaseGenotypeTable);
        PositionList goodPositions = new PositionListBuilder().addAll(goodPositionsFromNames).addAll(goodPositionsFromIndex).build();
        //Build the ones to remove, then compare
        PositionList badPositionsList = getPositionsFromSiteNames(mySiteNamesToRemove, myBaseGenotypeTable);
        PositionList siteConflicts = positionListsOverlap(goodPositionsFromNames, badPositionsList);
        if (siteConflicts.size() > 0) {
            StringBuilder intersects = new StringBuilder();
            for (Position p : siteConflicts) {
                intersects.append(p.getSNPID() + "\n");
            }
            throw new IllegalStateException("FilterGenotypeTableBuilder: build: Some sites manually specified both to keep and to remove:\n" + intersects.toString());
        }

        //Next compile a list of sites to remove. Will be inverted later, but lack of a remove() method makes it easier to accumulate the ones that fail
        PositionListBuilder badSiteBuilder = new PositionListBuilder();
        for(int s=0; s < myBaseGenotypeTable.numberOfSites(); s++){
            Position currentSite = myBaseGenotypeTable.positions().get(s);

            //Minor allele frequency
            double maf = myBaseGenotypeTable.minorAlleleFrequency(s);
            if(valueOutsideOfRange(maf, myMinFreqForSite, myMaxFreqForSite)){
                badSiteBuilder.add(currentSite);
                continue;
            }

            //Missingness
            //TODO: Include support for fractional values, not just counts
            int notMissing = myBaseGenotypeTable.totalNonMissingForSite(s);
            if(valueOutsideOfRange(notMissing, myMinCountForSite, myMaxCountForSite)){
                badSiteBuilder.add(currentSite);
                continue;
            }

            //Heterozygosity
            double hetCount = (double) myBaseGenotypeTable.heterozygousCount(s) / myBaseGenotypeTable.numberOfTaxa();
            if(valueOutsideOfRange(hetCount, myMinHeterozygousForSite, myMaxHeterozygousForSite)){
                badSiteBuilder.add(currentSite);
                continue;
            }
        }

        //Add in the ones specified to remove by name, then build
        badSiteBuilder.addAll(badPositionsList);
        PositionList badPositions = badSiteBuilder.build();

        //Go through and invert the list, returning all (good) sites that are EITHER specified to keep by name/index OR are not on the bad list
        PositionListBuilder goodSiteBuilder = new PositionListBuilder();
        for(Position p: myBaseGenotypeTable.positions()){
            if(goodPositionsFromNames.contains(p) || !badPositions.contains(p)){
                goodSiteBuilder.add(p);
            }
        }
        return goodSiteBuilder.build();
    }

    private boolean valueOutsideOfRange(double val, double min, double max){
        return (val < min || val > max);
    }

    public static TaxaList taxaListsOverlap(TaxaList a, TaxaList b) {
        TaxaListBuilder intersection = new TaxaListBuilder();
        for (Taxon t : a) {
            if (b.contains(t)) {
                intersection.add(t);
            }
        }
        return intersection.build();
    }

    //TODO: See if faster to change the list to something else first, eg, Set or Hash
    private static PositionList getPositionsFromSiteNames(List<String> names, GenotypeTable genos){
        PositionListBuilder posBuilder = new PositionListBuilder();
        for(int s=0; s<genos.numberOfSites(); s++){
            if(names.contains(genos.siteName(s))){
                posBuilder.add(genos.positions().get(s));
            }
        }
        return posBuilder.build();
    }

    private static PositionList getPositionsFromSiteIndex(int[] indices, GenotypeTable genos){
        PositionListBuilder posBuilder = new PositionListBuilder();
        for(int i: indices){
            posBuilder.add(genos.positions().get(i));
        }
        return posBuilder.build();
    }

    private static PositionList positionListsOverlap(PositionList a, PositionList b) {
        PositionListBuilder intersection = new PositionListBuilder();
        for (Position p : a) {
            if (b.contains(p)) {
                intersection.add(p);
            }
        }
        return intersection.build();
    }

}
