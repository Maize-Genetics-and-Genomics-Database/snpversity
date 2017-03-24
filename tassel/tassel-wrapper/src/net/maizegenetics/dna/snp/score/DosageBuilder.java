/*
 *  DosageBuilder
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.dna.snp.byte2d.FilterByte2D;
import net.maizegenetics.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class DosageBuilder {

    private Byte2DBuilder myBuilder;
    private final int myNumSites;

    private DosageBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        myBuilder = Byte2DBuilder.getInstance(numTaxa, numSites, SiteScore.SITE_SCORE_TYPE.Dosage, taxaList);
        myNumSites = numSites;
    }

    private DosageBuilder(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        myBuilder = Byte2DBuilder.getInstance(writer, numSites, SiteScore.SITE_SCORE_TYPE.Dosage, taxaList);
        myNumSites = numSites;
    }

    public static DosageBuilder getInstance(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        return new DosageBuilder(writer, numTaxa, numSites, taxaList);
    }

    public static DosageBuilder getInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new DosageBuilder(numTaxa, numSites, taxaList);
    }

    public static Dosage getFilteredInstance(Dosage base, FilterGenotypeTable filterGenotypeTable) {
        FilterByte2D resultStorage = Byte2DBuilder.getFilteredInstance(base.byteStorage(SiteScore.SITE_SCORE_TYPE.Dosage), filterGenotypeTable);
        return new Dosage(resultStorage);
    }

    public static Dosage getInstance(IHDF5Reader reader) {
        return new Dosage(Byte2DBuilder.getInstance(reader, SiteScore.SITE_SCORE_TYPE.Dosage));
    }

    public DosageBuilder addTaxon(int taxon, byte[] values) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("DosageBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        myBuilder.addTaxon(taxon, values);
        return this;
    }

    public Dosage build() {
        Byte2D input = myBuilder.build();
        myBuilder = null;
        return new Dosage(input);
    }

}
