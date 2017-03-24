package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;

import javax.activation.UnsupportedDataTypeException;
import javax.swing.*;
import java.awt.*;
import java.net.URL;

/**
 * Created by jgw87 on 6/5/14.
 * This plugin takes a Hapmap or VCF genotype file and sorts it according to TASSEL's conventions which rely on the
 * position, locus (chromosome), strand, and SNP name (to facilitate searching).
 */
public class SortGenotypeFilePlugin extends AbstractPlugin {

    private enum SupportedFileTypes {Hapmap, VCF}

    private PluginParameter<String> infile
            = new PluginParameter.Builder<>("inputFile", null, String.class)
            .required(true)
            .inFile()
            .guiName("Input file")
            .description("Input file")
            .build();
    private PluginParameter<String> outfile
            = new PluginParameter.Builder<>("outputFile", null, String.class)
            .required(true)
            .outFile()
            .guiName("Output file")
            .description("Output file")
            .build();
    private PluginParameter<SupportedFileTypes> fileType
            = new PluginParameter.Builder<>("fileType", null, SupportedFileTypes.class)
            .required(false)
            .guiName("File type")
            .description("Input/output file type (if not obvious from file name)")
            .build();


    public SortGenotypeFilePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet processData(DataSet input) {
        try {
            switch (fileType()) {
                case Hapmap:
                    GenotypeTable myHapmap = BuilderFromHapMap.getBuilder(inputFile(), null).sortPositions().build();
                    ExportUtils.writeToHapmap(myHapmap, outputFile());
                    break;
                case VCF:
                    GenotypeTable myVCF = BuilderFromVCF.getBuilder(inputFile()).keepDepth().buildAndSortInMemory();
                    ExportUtils.writeToVCF(myVCF, outputFile(), true);
                    break;
                default:
                    throw new UnsupportedDataTypeException("SortGenotypeFilePlugin: Identified data type does not conform to known types (Hapmap, VCF)");
            }
        } catch (UnsupportedDataTypeException e) {
            e.printStackTrace();
        }

        return null;
    }

    @Override
    public String pluginDescription() {
        return "This plugin takes a Hapmap or VCF genotype file and sorts it according to TASSEL's conventions, " +
                "which rely on the position, locus (chromosome), strand, and SNP name (to facilitate searching).";
    }

    @Override
    protected void postProcessParameters() {
        //If file type not provided, try to guess from file name
        if (fileType() == null) {
            if (inputFile().toLowerCase().endsWith(".hmp.txt") || inputFile().toLowerCase().endsWith(".hmp.txt.gz")) {
                fileType(SupportedFileTypes.Hapmap);
            } else if (inputFile().toLowerCase().endsWith(".vcf") || inputFile().toLowerCase().endsWith(".vcf.gz")) {
                fileType(SupportedFileTypes.VCF);
            } else {
                throw new UnsupportedOperationException("Unable to guess file type from input file name. Please rename to end in .hmp.txt, .hmp.txt.gz, .vcf, or .vcf.gz");
            }
        }
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SortGenotypeFilePlugin.class.getResource("/net/maizegenetics/analysis/images/sort.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Sort Genotype File";
    }

    @Override
    public String getToolTipText() {
        return "Sort Genotype File";
    }


    public static void main(String[] args) {
        GeneratePluginCode.generate(SortGenotypeFilePlugin.class);
    }

    //Code below this point was automatically generated with the above main() method. Do not manually alter

    /**
     * Input file
     *
     * @return Input file
     */
    public String inputFile() {
        return infile.value();
    }

    /**
     * Set Input file. Input file
     *
     * @param value Input file
     * @return this plugin
     */
    public SortGenotypeFilePlugin inputFile(String value) {
        infile = new PluginParameter<>(infile, value);
        return this;
    }

    /**
     * Output file
     *
     * @return Output file
     */
    public String outputFile() {
        return outfile.value();
    }

    /**
     * Set Output file. Output file
     *
     * @param value Output file
     * @return this plugin
     */
    public SortGenotypeFilePlugin outputFile(String value) {
        outfile = new PluginParameter<>(outfile, value);
        return this;
    }

    /**
     * Input/output file type (if not obvious from file name)
     *
     * @return File type
     */
    public SupportedFileTypes fileType() {
        return fileType.value();
    }

    /**
     * Set File type. Input/output file type (if not obvious
     * from file name)
     *
     * @param value File type
     * @return this plugin
     */
    public SortGenotypeFilePlugin fileType(SupportedFileTypes value) {
        fileType = new PluginParameter<>(fileType, value);
        return this;
    }
}
