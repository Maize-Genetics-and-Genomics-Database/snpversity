/*
 *  AlleleDepthBuilder
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.dna.snp.byte2d.FilterByte2D;
import net.maizegenetics.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleDepthBuilder {

    private final Map<SiteScore.SITE_SCORE_TYPE, Byte2DBuilder> myBuilders = new LinkedHashMap<>();
    private final int myNumSites;

    private AlleleDepthBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        for (int i = 0; i < AlleleDepth.ALLELE_DEPTH_TYPES.length; i++) {
            myBuilders.put(AlleleDepth.ALLELE_DEPTH_TYPES[i], Byte2DBuilder.getInstance(numTaxa, numSites, AlleleDepth.ALLELE_DEPTH_TYPES[i], taxaList));
        }
        myNumSites = numSites;
    }

    private AlleleDepthBuilder(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        for (int i = 0; i < AlleleDepth.ALLELE_DEPTH_TYPES.length; i++) {
            myBuilders.put(AlleleDepth.ALLELE_DEPTH_TYPES[i], Byte2DBuilder.getInstance(writer, numSites, AlleleDepth.ALLELE_DEPTH_TYPES[i], taxaList));
        }
        myNumSites = numSites;
    }

    /**
     * This returns an AlleleDepthBuilder where depths are stored in a HDF5
     * file.
     *
     * @param writer writer
     * @param numTaxa number of taxa
     * @param numSites number of sites
     * @param taxaList taxa list
     *
     * @return AlleleDepthBuilder
     */
    public static AlleleDepthBuilder getInstance(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleDepthBuilder(writer, numTaxa, numSites, taxaList);
    }

    /**
     * This returns an AlleleDepthBuilder where depths are stored in memory.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites
     * @param taxaList taxa list
     *
     * @return AlleleDepthBuilder
     */
    public static AlleleDepthBuilder getInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleDepthBuilder(numTaxa, numSites, taxaList);
    }

    /**
     * This creates a filtered AlleleDepth.
     *
     * @param base original AlleleDepth
     * @param filterGenotypeTable filter
     *
     * @return filtered AlleleDepth
     */
    public static AlleleDepth getFilteredInstance(AlleleDepth base, FilterGenotypeTable filterGenotypeTable) {
        Collection<Byte2D> storage = base.byteStorage();
        FilterByte2D[] resultStorage = new FilterByte2D[storage.size()];
        int count = 0;
        for (Byte2D current : storage) {
            resultStorage[count] = Byte2DBuilder.getFilteredInstance(current, filterGenotypeTable);
        }
        return new AlleleDepth(resultStorage);
    }

    /**
     * This creates an AlleleDepth instance from an existing HDF5 file.
     *
     * @param reader reader
     *
     * @return AlleleDepth
     */
    public static AlleleDepth getInstance(IHDF5Reader reader) {
        int numAlleles = AlleleDepth.ALLELE_DEPTH_TYPES.length;
        Byte2D[] input = new Byte2D[numAlleles];
        int count = 0;
        for (SiteScore.SITE_SCORE_TYPE current : AlleleDepth.ALLELE_DEPTH_TYPES) {
            input[count++] = Byte2DBuilder.getInstance(reader, current);
        }
        return new AlleleDepth(input);
    }

    public AlleleDepthBuilder addTaxon(int taxon, int[] values, SiteScore.SITE_SCORE_TYPE type) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("AlleleDepthBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        byte[] result = AlleleDepthUtil.depthIntToByte(values);
        myBuilders.get(type).addTaxon(taxon, result);
        return this;
    }

    public AlleleDepthBuilder addTaxon(int taxon, byte[][] values) {
        if (AlleleDepth.ALLELE_DEPTH_TYPES.length != values.length) {
            throw new IllegalArgumentException("AlleleDepthBuilder: addTaxon: number of alleles: " + values.length + " doesn't equals: " + AlleleDepth.ALLELE_DEPTH_TYPES.length);
        }
        if (myNumSites != values[0].length) {
            throw new IllegalArgumentException("AlleleDepthBuilder: addTaxon: number of values: " + values[0].length + " doesn't equal number of sites: " + myNumSites);
        }
        int count = 0;
        for (SiteScore.SITE_SCORE_TYPE current : AlleleDepth.ALLELE_DEPTH_TYPES) {
            myBuilders.get(current).addTaxon(taxon, values[count++]);
        }
        return this;
    }

    public AlleleDepth build() {
        Byte2D[] input = new Byte2D[myBuilders.size()];
        int count = 0;
        for (Byte2DBuilder builder : myBuilders.values()) {
            input[count++] = builder.build();
        }
        myBuilders.clear();
        return new AlleleDepth(input);
    }

}
