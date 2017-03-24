/*
 * GenotypeTableMaskGeneticDistance
 */
package net.maizegenetics.gui;

import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.dna.snp.GenotypeTable;

import java.awt.*;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeTableMaskGeneticDistance extends AbstractGenotypeTableMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private Map<Integer, Byte> myCache = new LinkedHashMap<Integer, Byte>() {
        protected boolean removeEldestEntry(Map.Entry<Integer, Byte> eldest) {
            return size() > 100;
        }
    };
    private final int myTaxonReference;
    private final GenotypeTable myAlignment;

    private GenotypeTableMaskGeneticDistance(GenotypeTable align, int taxonReference, String name, Color color) {
        super(align, name, color, GenotypeTableMask.MaskType.reference);
        myTaxonReference = taxonReference;
        myAlignment = align;
    }

    public static GenotypeTableMaskGeneticDistance getInstanceCompareReference(GenotypeTable align, Taxon id) {
        int index = align.taxa().indexOf(id);
        if (index < 0) {
            throw new IllegalArgumentException("GenotypeTableMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static GenotypeTableMaskGeneticDistance getInstanceCompareReference(GenotypeTable align, String id) {
        int index = align.taxa().indexOf(id);
        if (index < 0) {
            throw new IllegalArgumentException("GenotypeTableMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static GenotypeTableMaskGeneticDistance getInstanceCompareReference(GenotypeTable align, int index) {
        if ((index < 0) || (index >= align.numberOfTaxa())) {
            throw new IllegalArgumentException("GenotypeTableMaskGeneticDistance: getInstanceCompareReference: unknown index: " + index);
        }
        String name = align.taxaName(index) + " Genetic Distance";
        return new GenotypeTableMaskGeneticDistance(align, index, name, null);
    }

    @Override
    public byte getMask(int taxon, int site) {

        Byte result = myCache.get(taxon);
        if (result != null) {
            return result;
        }

        result = (byte) (IBSDistanceMatrix.computeHetBitDistances(myAlignment, taxon, myTaxonReference)[0] * 255.0);
        myCache.put(taxon, result);
        return result;

    }
}
