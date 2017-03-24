package net.maizegenetics.analysis.distance;

import com.google.common.collect.Range;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.distance.DistanceMatrix;

import javax.swing.*;

import java.net.URL;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * @author Terry Casstevens
 * @author Zhiwu Zhang
 * @author Peter Bradbury
 *
 */
public class KinshipPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(KinshipPlugin.class);

    public static enum KINSHIP_METHOD {

        Centered_IBS,
        Normalized_IBS,
        Dominance_Centered_IBS,
        Dominance_Normalized_IBS
    };

    private final GenotypeTable.GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GenotypeTable.GENOTYPE_TABLE_COMPONENT[]{
        GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability};

    private PluginParameter<KINSHIP_METHOD> myMethod = new PluginParameter.Builder<>("method", KINSHIP_METHOD.Centered_IBS, KINSHIP_METHOD.class)
            .guiName("Kinship method")
            .range(KINSHIP_METHOD.values())
            .description("The Centered_IBS (Endelman - previously Scaled_IBS) method produces a kinship matrix that is scaled to give a reasonable estimate of additive "
                    + "genetic variance. Uses algorithm http://www.g3journal.org/content/2/11/1405.full.pdf Equation-13. "
                    + "The Normalized_IBS (Previously GCTA) uses the algorithm published here: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/pdf/main.pdf.")
            .build();

    private PluginParameter<Integer> myMaxAlleles = new PluginParameter.Builder<>("maxAlleles", 6, Integer.class)
            .description("")
            .range(Range.closed(2, 6))
            .dependentOnParameter(myMethod, new Object[]{KINSHIP_METHOD.Centered_IBS, KINSHIP_METHOD.Dominance_Centered_IBS})
            .build();

    public KinshipPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if ((alignInList == null) || (alignInList.isEmpty())) {
            throw new IllegalArgumentException("KinshipPlugin: Nothing selected. Please select a genotype.");
        }
    }

    @Override
    public void setParameters(String[] args) {
        try {
            for (int i = 0; i < args.length; i++) {
                if (args[i].equals("-method")) {
                    if (args[i + 1].equalsIgnoreCase("GCTA")) {
                        args[i + 1] = KINSHIP_METHOD.Normalized_IBS.name();
                        myLogger.warn("setParameters: Notice GCTA has been changed to Normalized_IBS");
                    } else if (args[i + 1].equalsIgnoreCase("Scaled_IBS")) {
                        args[i + 1] = KINSHIP_METHOD.Centered_IBS.name();
                        myLogger.warn("setParameters: Notice Scaled_IBS has been changed to Centered_IBS");
                    } else if (args[i + 1].equalsIgnoreCase("Dominance")) {
                        args[i + 1] = KINSHIP_METHOD.Dominance_Centered_IBS.name();
                        myLogger.warn("setParameters: Notice Dominance has been changed to Dominance_Centered_IBS");
                    }
                }
                break;
            }
        } catch (Exception e) {
            // do nothing
            myLogger.debug(e.getMessage(), e);
        }
        super.setParameters(args);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

        List<Datum> result = new ArrayList<>();
        Iterator<Datum> itr = alignInList.iterator();
        while (itr.hasNext()) {

            Datum current = itr.next();
            String datasetName = current.getName();
            DistanceMatrix kin = null;

            if (current.getData() instanceof GenotypeTable) {
                GenotypeTable myGenotype = (GenotypeTable) current.getData();
                if (kinshipMethod() == KINSHIP_METHOD.Centered_IBS) {
                    kin = EndelmanDistanceMatrix.getInstance(myGenotype, maxAlleles(), this);
                } else if (kinshipMethod() == KINSHIP_METHOD.Normalized_IBS) {
                    kin = GCTADistanceMatrix.getInstance(myGenotype, this);
                } else if (kinshipMethod() == KINSHIP_METHOD.Dominance_Centered_IBS) {
                    kin = DominanceRelationshipMatrix.getInstance(myGenotype, maxAlleles(), this);
                } else if (kinshipMethod() == KINSHIP_METHOD.Dominance_Normalized_IBS) {
                    throw new UnsupportedOperationException("Method Dominance_Normalized_IBS hasn't been implemented yet.");
                } else {
                    throw new IllegalArgumentException("Unknown method to calculate kinship: " + kinshipMethod());
                }
            } else {
                throw new IllegalArgumentException("Invalid selection. Can't create kinship matrix from: " + datasetName);
            }

            if (kin != null) {
                //add kin to datatree;
                Datum ds = new Datum(kinshipMethod() + "_" + datasetName, kin, kinshipMethod() + " matrix created from " + datasetName);
                result.add(ds);
            }

        }

        return new DataSet(result, this);

    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = KinshipPlugin.class.getResource("/net/maizegenetics/analysis/images/Kin.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Kinship";
    }

    @Override
    public String getToolTipText() {
        return "Calculate kinship from marker data";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(KinshipPlugin.class);
    // }
    /**
     * Convenience method to run plugin with one return object.
     */
    public DistanceMatrix runPlugin(DataSet input) {
        return (DistanceMatrix) performFunction(input).getData(0).getData();
    }

    /**
     * The scaled_IBS method produces a kinship matrix that is scaled to give a
     * reasonable estimate of additive genetic variance. The pairwise_IBS
     * method, which is the method used by TASSEL ver.4, may result in an
     * inflated estimate of genetic variance. Either will do a good job of
     * controlling population structure in MLM. The pedigree method is used to
     * calculate a kinship matrix from a pedigree information.
     *
     * @return Kinship method
     */
    public KINSHIP_METHOD kinshipMethod() {
        return myMethod.value();
    }

    /**
     * Set Kinship method. The scaled_IBS method produces a kinship matrix that
     * is scaled to give a reasonable estimate of additive genetic variance. The
     * pairwise_IBS method, which is the method used by TASSEL ver.4, may result
     * in an inflated estimate of genetic variance. Either will do a good job of
     * controlling population structure in MLM. The pedigree method is used to
     * calculate a kinship matrix from a pedigree information.
     *
     * @param value Kinship method
     *
     * @return this plugin
     */
    public KinshipPlugin kinshipMethod(KINSHIP_METHOD value) {
        myMethod = new PluginParameter<>(myMethod, value);
        return this;
    }

    /**
     * Max Alleles
     *
     * @return Max Alleles
     */
    public Integer maxAlleles() {
        return myMaxAlleles.value();
    }

    /**
     * Set Max Alleles. Max Alleles
     *
     * @param value Max Alleles
     *
     * @return this plugin
     */
    public KinshipPlugin maxAlleles(Integer value) {
        myMaxAlleles = new PluginParameter<>(myMaxAlleles, value);
        return this;
    }

}
