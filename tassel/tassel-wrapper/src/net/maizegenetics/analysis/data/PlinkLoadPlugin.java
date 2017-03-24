/*
 * PlinkLoadPlugin.java
 *
 * Created on December 18, 2009
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

import javax.swing.*;

import java.awt.Frame;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class PlinkLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(PlinkLoadPlugin.class);

    private PluginParameter<String> myPedFile = new PluginParameter.Builder<>("pedFile", null, String.class).required(true).inFile()
            .description("Ped File").build();
    private PluginParameter<String> myMapFile = new PluginParameter.Builder<>("mapFile", null, String.class).required(true).inFile()
            .description("Map File").build();

    /**
     * Creates a new instance of PlinkLoadPlugin
     */
    public PlinkLoadPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet processData(DataSet input) {
        return loadFile(pedFile(), mapFile());
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(PlinkLoadPlugin.class);
    // }
    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Ped File
     *
     * @return Ped File
     */
    public String pedFile() {
        return myPedFile.value();
    }

    /**
     * Set Ped File. Ped File
     *
     * @param value Ped File
     *
     * @return this plugin
     */
    public PlinkLoadPlugin pedFile(String value) {
        myPedFile = new PluginParameter<>(myPedFile, value);
        return this;
    }

    /**
     * Map File
     *
     * @return Map File
     */
    public String mapFile() {
        return myMapFile.value();
    }

    /**
     * Set Map File. Map File
     *
     * @param value Map File
     *
     * @return this plugin
     */
    public PlinkLoadPlugin mapFile(String value) {
        myMapFile = new PluginParameter<>(myMapFile, value);
        return this;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Load Plink";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load Plink Files";
    }

    public DataSet loadFile(String thePedFile, String theMapFile) {

        GenotypeTable result = ImportUtils.readFromPLink(thePedFile, theMapFile, this);
        Datum td = new Datum(Utils.getFilename(thePedFile, FileLoadPlugin.FILE_EXT_PLINK_PED), result, null);
        DataSet tds = new DataSet(td, this);
        fireDataSetReturned(new PluginEvent(tds, PlinkLoadPlugin.class));

        return tds;

    }

}
