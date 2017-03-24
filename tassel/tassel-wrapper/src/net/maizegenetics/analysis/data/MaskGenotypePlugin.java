/*
 *  MaskGenotypePlugin
 * 
 *  Created on May 8, 2015
 */
package net.maizegenetics.analysis.data;

import com.google.common.collect.Range;
import java.awt.Frame;
import javax.swing.*;
import java.util.List;
import java.util.Random;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.MaskedGenotypes;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;

/**
 * The purpose of this class is to mask (make UNKNOWN) a portion of the
 * genotypes of the given GenotypeTable. The resulting GenotypeTable will be the
 * same except some genotypes will have been changed to Unknown.
 *
 * @author Terry Casstevens
 */
public class MaskGenotypePlugin extends AbstractPlugin {

    private PluginParameter<Double> myPercentageMasked
            = new PluginParameter.Builder<>("percentageMasked", 0.0, Double.class)
            .range(Range.closed(0.0, 1.0))
            .description("Percentage of genotypes to mask.")
            .build();

    public MaskGenotypePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> data = input.getDataOfType(GenotypeTable.class);
        if (data.size() != 1) {
            throw new IllegalArgumentException("MaskGenotypePlugin: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum inputDatum = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable original = (GenotypeTable) inputDatum.getData();

        GenotypeCallTable origCalls = original.genotypeMatrix();

        GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getInstanceCopy(origCalls);
        MaskedGenotypes mask = new MaskedGenotypes(builder.getTaxaCount(), builder.getSiteCount());

        percentageMask(builder, origCalls, mask);

        GenotypeTable result = GenotypeTableBuilder.getInstance(original, builder.build());

        Datum genotype = new Datum(inputDatum.getName() + "_Masked", result, null);
        Datum maskDatum = new Datum(inputDatum.getName() + "_Mask", mask, null);

        return new DataSet(new Datum[]{genotype, maskDatum}, this);

    }

    private void percentageMask(GenotypeCallTableBuilder builder, GenotypeCallTable origCalls, MaskedGenotypes mask) {
        if (percentageMasked() == 0.0) {
            return;
        }

        Random random = new Random(1);
        int numTaxa = builder.getTaxaCount();
        int numSites = builder.getSiteCount();
        long totalGenotypes = numTaxa * numSites;
        long numKnown = origCalls.stream().filter(b -> b != GenotypeTable.UNKNOWN_DIPLOID_ALLELE).count();
        double percentKnown = numKnown / totalGenotypes;
        if (percentKnown <= percentageMasked()) {
            throw new IllegalArgumentException("MaskGenotypePlugin: percentMask: Trying to mask: " + percentageMasked() + "%.  But this genotype only has: " + percentKnown + "% Known.");
        }
        long numberMasked = 0;
        long numberToMask = (long) (totalGenotypes / percentageMasked());
        long step = (long) (1.0 / percentageMasked() * 2.0);
        long site = 0;
        while (true) {
            for (int t = 0; t < numTaxa;) {
                site += (long) ((double) step * random.nextDouble());
                if (site >= numSites) {
                    site %= numSites;
                    t += site / numSites;
                }
                if ((!mask.get(t, (int) site))
                        && (origCalls.genotype(t, (int) site) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)) {
                    builder.setBase(t, (int) site, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                    mask.set(t, (int) site);
                    numberMasked++;
                    if (numberMasked >= numberToMask) {
                        return;
                    }
                }
            }
        }
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Percentage of genotypes to mask.
     *
     * @return Percentage Masked
     */
    public Double percentageMasked() {
        return myPercentageMasked.value();
    }

    /**
     * Set Percentage Masked. Percentage of genotypes to mask.
     *
     * @param value Percentage Masked
     *
     * @return this plugin
     */
    public MaskGenotypePlugin percentageMasked(Double value) {
        myPercentageMasked = new PluginParameter<>(myPercentageMasked, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Mask Genotype";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Mask Genotype";
    }

}
