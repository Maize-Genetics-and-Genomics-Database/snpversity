package net.maizegenetics.tassel;

import java.awt.Frame;
import javax.swing.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

import net.maizegenetics.prefs.TasselPrefs;

/**
 * @author Terry Casstevens
 */
public class PreferencesDialog extends AbstractPlugin {

    private PluginParameter<Boolean> myRetainRareAlleles = new PluginParameter.Builder<>("retainRareAlleles", TasselPrefs.ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT, Boolean.class)
            .description("True if rare alleles should be retained.  This has no effect on Nucleotide Data as all alleles will be retained regardless.")
            .build();

    public PreferencesDialog(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        setParameter(myRetainRareAlleles, TasselPrefs.getAlignmentRetainRareAlleles());
    }

    @Override
    public DataSet processData(DataSet input) {
        TasselPrefs.putAlignmentRetainRareAlleles(retainRareAlleles());
        return null;
    }

    /**
     * Retain Rare Alleles
     *
     * @return Retain Rare Alleles
     */
    public Boolean retainRareAlleles() {
        return myRetainRareAlleles.value();
    }

    /**
     * Set Retain Rare Alleles. Retain Rare Alleles
     *
     * @param value Retain Rare Alleles
     *
     * @return this plugin
     */
    public PreferencesDialog retainRareAlleles(Boolean value) {
        myRetainRareAlleles = new PluginParameter<>(myRetainRareAlleles, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Preferences";
    }

    @Override
    public String getToolTipText() {
        return "Preferences";
    }
}
