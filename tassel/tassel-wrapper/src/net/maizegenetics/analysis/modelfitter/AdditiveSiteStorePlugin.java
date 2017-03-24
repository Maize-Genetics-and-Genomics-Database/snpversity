package net.maizegenetics.analysis.modelfitter;

import java.awt.Frame;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.analysis.modelfitter.AdditiveSite.CRITERION;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;

public class AdditiveSiteStorePlugin extends AbstractPlugin {
    private static Logger myLogger = Logger.getLogger(AdditiveSiteStorePlugin.class);
    PluginParameter<String> outFilename =
            new PluginParameter.Builder<>("filename", null, String.class)
                    .guiName("Output Filename")
                    .description("The name of file to which the additive sites will be written.")
                    .required(true)
                    .build();
    PluginParameter<Boolean> isGenotype =
            new PluginParameter.Builder<>("genotype", true, Boolean.class)
                    .guiName("Is Genotype?")
                    .description("Is the input nucleotide genotypes? If false, then ReferenceProbability will be used.")
                    .build();

    public AdditiveSiteStorePlugin() {
        super(null, false);
    }

    public AdditiveSiteStorePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> datumList = input.getDataOfType(GenotypeTable.class);
        if (datumList.size() != 0)
            throw new IllegalArgumentException("Exactly one genotype dataset will is required as input.");
    }

    @Override
    public DataSet processData(DataSet input) {
        GenotypeTable myGenotype =
                (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
        int numberOfSites = myGenotype.numberOfSites();

        List<AdditiveSite> siteList = new ArrayList<>(numberOfSites);
        long start = System.nanoTime();
        if (isGenotype.value()) {
            for (int s = 0; s < numberOfSites; s++) {
                siteList.add(new GenotypeAdditiveSite(s, myGenotype.chromosomeName(s), myGenotype.chromosomalPosition(s), myGenotype.siteName(s), CRITERION.pval, myGenotype.genotypeAllTaxa(s), myGenotype.majorAllele(s), myGenotype.majorAlleleFrequency(s)));
            }
        } else {
            int ntaxa = myGenotype.numberOfTaxa();
            for (int s = 0; s < numberOfSites; s++) {
                float[] refprobs = new float[ntaxa];
                for (int t = 0; t < ntaxa; t++)
                    refprobs[t] = myGenotype.referenceProbability(t, s);
                siteList.add(new RefProbAdditiveSite(s, myGenotype.chromosomeName(s), myGenotype.chromosomalPosition(s), myGenotype.siteName(s), CRITERION.pval, refprobs));
            }

        }

        myLogger.debug(String.format("site list created with %d sites in %d ms.", siteList.size(), (System.nanoTime() - start) / 1000000));

        start = System.nanoTime();
        String outfile = outFilename.value();
        ObjectOutputStream out;
        try {
            out = new ObjectOutputStream(new FileOutputStream(outfile));

            Consumer<AdditiveSite> writeSite = s -> {
                try {
                    out.writeObject(s);
                } catch (IOException ioe) {
                    throw new RuntimeException("Error writing to SiteList object store.", ioe);
                } 
            };

            out.writeObject(new Integer(siteList.size()));
            out.writeObject(myGenotype.taxa());
            siteList.stream().forEach(writeSite);
            out.flush();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        myLogger.debug(String.format("Sites written to outfile in %d ms.\n", (System.nanoTime() - start) / 1000000));

        return null;
    }

    @Override
    public String pluginDescription() {
        return "This plugin serializes a genotype as a list of additive sites. The objects stored are an Integer equal to the number of sites, "
                + "the Genotype TaxaList, then all of the additive sites";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Write additive site file";
    }

    @Override
    public String getToolTipText() {
        return "Stores genotypes coded as additive sites";
    }

}
