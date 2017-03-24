/*
 * ExportPlugin.java
 *
 * Created on December 18, 2009
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.FilterList;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.FilterJSONUtils;
import net.maizegenetics.dna.snp.io.SiteScoresIO;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListTableReport;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListTableReport;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.WriteDistanceMatrix;
import net.maizegenetics.taxa.tree.SimpleTree;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.util.*;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;

/**
 *
 * @author Terry Casstevens
 */
public class ExportPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ExportPlugin.class);
    private FileLoadPlugin.TasselFileType myFileType = FileLoadPlugin.TasselFileType.Hapmap;
    private String mySaveFile = null;
    private boolean myKeepDepth = true;
    private boolean myIncludeTaxaAnnotations = TasselPrefs.EXPORT_PLUGIN_INCLUDE_TAXA_ANNOTATIONS_DEFAULT;
    private SiteScore.SITE_SCORE_TYPE mySiteScoreType = null;
    private final JFileChooser myFileChooserSave;

    /**
     * Creates a new instance of ExportPlugin
     */
    public ExportPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
        if (isInteractive) {
            myFileChooserSave = new JFileChooser(TasselPrefs.getSaveDir());
            myFileChooserSave.setMultiSelectionEnabled(false);
        } else {
            myFileChooserSave = null;
        }
    }

    public DataSet performFunction(DataSet input) {

        try {

            if (input.getSize() != 1) {
                String message = "Please select one and only one item.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error("performFunction: " + message);
                }
                return null;
            }

            String filename = mySaveFile;
            try {
                Object data = input.getData(0).getData();
                if (data instanceof GenotypeTable) {
                    filename = performFunctionForAlignment((GenotypeTable) data);
                } else if (data instanceof Phenotype) {
                    filename = performFunctionForPhenotype((Phenotype) data);
                } else if (data instanceof FilterList) {
                    filename = performFunctionForFilter((FilterList) data);
                } else if (data instanceof DistanceMatrix) {
                    filename = performFunctionForDistanceMatrix((DistanceMatrix) data);
                } else if (data instanceof TaxaList) {
                    filename = performFunctionForTaxaList((TaxaList) data);
                } else if (data instanceof TaxaListTableReport) {
                    filename = performFunctionForTaxaList(((TaxaListTableReport) data).getTaxaList());
                } else if (data instanceof PositionList) {
                    filename = performFunctionForPositionList((PositionList) data);
                } else if (data instanceof PositionListTableReport) {
                    filename = performFunctionForPositionList(((PositionListTableReport) data).getPositionList());
                } else if (data instanceof TableReport) {
                    filename = performFunctionForTableReport((TableReport) data);
                } else if (data instanceof SimpleTree) {
                    filename = performFunctionForSimpleTree((SimpleTree) data);
                } else {
                    String message = "Don't know how to export data type: " + data.getClass().getName();
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), message);
                    } else {
                        myLogger.error("performFunction: " + message);
                    }
                    return null;
                }
            } catch (Exception e) {
                e.printStackTrace();
                StringBuilder builder = new StringBuilder();
                builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
                String str = builder.toString();
                if (isInteractive()) {
                    DialogUtils.showError(str, getParentFrame());
                } else {
                    myLogger.error(str);
                }

                return null;
            }

            if (filename != null) {
                myLogger.info("performFunction: wrote dataset: " + input.getData(0).getName() + " to file: " + filename);
                return new DataSet(new Datum("Filename", filename, null), this);
            } else {
                return null;
            }

        } finally {
            fireProgress(100);
        }

    }

    public String performFunctionForDistanceMatrix(DistanceMatrix input) {

        if (isInteractive()) {
            ExportSquareMatrixDialog theDialog = new ExportSquareMatrixDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }
            myFileType = theDialog.getTasselFileType();
            theDialog.dispose();
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        try {

            if (myFileType != FileLoadPlugin.TasselFileType.SqrMatrixRaw && myFileType != FileLoadPlugin.TasselFileType.SqrMatrixBin) {
                String filename = Utils.addSuffixIfNeeded(mySaveFile, ".txt", new String[]{".txt", ".txt.gz"});
                WriteDistanceMatrix.saveDelimitedDistanceMatrix(input, filename);
                return filename;
            } else if (myFileType == FileLoadPlugin.TasselFileType.SqrMatrixRaw) {
                File taxaListFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".grm.id"));
                File matrixFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".grm.raw"));
                WriteDistanceMatrix.saveRawMultiBlupMatrix(input, taxaListFile, matrixFile);
                return matrixFile.getCanonicalPath();
            } else if (myFileType == FileLoadPlugin.TasselFileType.SqrMatrixBin) {
                File taxaListFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".grm.id"));
                File countsFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".grm.N.bin"));
                File kinshipFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".grm.bin"));
                WriteDistanceMatrix.saveBinMultiBlupMatrix(input, taxaListFile, kinshipFile, countsFile);
                return kinshipFile.getCanonicalPath();
            } else {
                throw new IllegalArgumentException("ExportPlugin: performFunctionForDistanceMatrix: Unknown file type: " + myFileType);
            }

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("ExportPlugin: performFunctionForDistanceMatrix: Problem writing file: " + mySaveFile);
        }

    }

    public String performFunctionForTableReport(TableReport input) {
        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        try {
            File theFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".txt"));
            TableReportUtils.saveDelimitedTableReport(input, "\t", theFile);
            return theFile.getCanonicalPath();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("ExportPlugin: performFunctionForTableReport: Problem writing file: " + mySaveFile);
        }

    }

    public String performFunctionForFilter(FilterList filter) {

        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        try {
            return FilterJSONUtils.exportFilterToJSON(filter, mySaveFile);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ExportPlugin: performFunctionForFilter: Problem writing file: " + mySaveFile + "\n" + e.getMessage());
        }

    }

    public String performFunctionForPhenotype(Phenotype input) {

        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String filename = "";
        try {
            filename = Utils.addSuffixIfNeeded(mySaveFile, ".txt");
            PhenotypeUtils.write(input, filename);
            return new File(filename).getCanonicalPath();
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ExportPlugin: performFunctionForPhenotype: Problem writing file: " + filename);
        }

    }

    public String performFunctionForAlignment(GenotypeTable inputAlignment) {

        java.util.List<GenotypeTable.GENOTYPE_TABLE_COMPONENT> components = new ArrayList<>();
        if (inputAlignment.hasGenotype()) {
            components.add(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype);
        }
        if (inputAlignment.hasReferenceProbablity()) {
            components.add(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability);
        }
        if (inputAlignment.hasAlleleProbabilities()) {
            components.add(GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability);
        }
        if (inputAlignment.hasDepth()) {
            components.add(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Depth);
        }
        if (inputAlignment.hasDosage()) {
            components.add(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Dosage);
        }

        if (isInteractive()) {

            myFileType = null;
            mySiteScoreType = null;
            if ((components.contains(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability)) && (components.size() == 1)) {
                mySiteScoreType = SiteScore.SITE_SCORE_TYPE.ReferenceProbablity;
            } else {
                ExportPluginDialog theDialog = new ExportPluginDialog(components);
                theDialog.setLocationRelativeTo(getParentFrame());
                theDialog.setVisible(true);
                if (theDialog.isCancel()) {
                    return null;
                }
                myFileType = theDialog.getTasselFileType();
                if (myFileType == null) {
                    mySiteScoreType = theDialog.getSiteScoreType();
                }
                myKeepDepth = theDialog.keepDepth();

                theDialog.dispose();
            }

            setSaveFile(getFileByChooser());
        } else if ((components.contains(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability)) && (components.size() == 1)) {
            mySiteScoreType = SiteScore.SITE_SCORE_TYPE.ReferenceProbablity;
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String resultFile = mySaveFile;

        if (mySiteScoreType == SiteScore.SITE_SCORE_TYPE.ReferenceProbablity) {
            resultFile = SiteScoresIO.writeReferenceProbability(inputAlignment, resultFile);
        } else if ((myFileType == FileLoadPlugin.TasselFileType.Hapmap) || (myFileType == FileLoadPlugin.TasselFileType.HapmapDiploid)) {
            boolean isDiploid = false;
            if (isInteractive()) {
                HapmapOptionDialog diploidDialog = new HapmapOptionDialog();
                diploidDialog.setLocationRelativeTo(getParentFrame());
                diploidDialog.setVisible(true);
                if (diploidDialog.isCancel()) {
                    return null;
                }
                isDiploid = diploidDialog.getDiploid();
                myIncludeTaxaAnnotations = diploidDialog.includeTaxaAnnotations();
            } else {
                if (myFileType == FileLoadPlugin.TasselFileType.Hapmap) {
                    isDiploid = false;
                } else if (myFileType == FileLoadPlugin.TasselFileType.HapmapDiploid) {
                    isDiploid = true;
                }
            }
            resultFile = ExportUtils.writeToHapmap(inputAlignment, isDiploid, mySaveFile, '\t', myIncludeTaxaAnnotations, this);
        } else if (myFileType == FileLoadPlugin.TasselFileType.Plink) {
            resultFile = ExportUtils.writeToPlink(inputAlignment, mySaveFile, '\t');
        } else if (myFileType == FileLoadPlugin.TasselFileType.Phylip_Seq) {
            PrintWriter out = null;
            try {
                resultFile = Utils.addSuffixIfNeeded(mySaveFile, ".phy");
                out = new PrintWriter(new FileWriter(resultFile));
                ExportUtils.printSequential(inputAlignment, out);
            } catch (Exception e) {
                throw new IllegalStateException("ExportPlugin: performFunction: Problem writing file: " + mySaveFile);
            } finally {
                out.flush();
                out.close();
            }
        } else if (myFileType == FileLoadPlugin.TasselFileType.Phylip_Inter) {
            PrintWriter out = null;
            try {
                resultFile = Utils.addSuffixIfNeeded(mySaveFile, ".phy");
                out = new PrintWriter(new FileWriter(resultFile));
                ExportUtils.printInterleaved(inputAlignment, out);
            } catch (Exception e) {
                throw new IllegalStateException("ExportPlugin: performFunction: Problem writing file: " + mySaveFile);
            } finally {
                out.flush();
                out.close();
            }
        } else if (myFileType == FileLoadPlugin.TasselFileType.Table) {
            resultFile = ExportUtils.saveDelimitedAlignment(inputAlignment, "\t", mySaveFile);
        } else if (myFileType == FileLoadPlugin.TasselFileType.Serial) {
            resultFile = ExportUtils.writeAlignmentToSerialGZ(inputAlignment, mySaveFile);
        } else if (myFileType == FileLoadPlugin.TasselFileType.HDF5) {
            resultFile = ExportUtils.writeGenotypeHDF5(inputAlignment, mySaveFile, myKeepDepth);
        } else if (myFileType == FileLoadPlugin.TasselFileType.VCF) {
            resultFile = ExportUtils.writeToVCF(inputAlignment, mySaveFile, myKeepDepth);
        } else {
            throw new IllegalStateException("ExportPlugin: performFunction: Unknown Alignment File Format: " + myFileType);
        }

        return resultFile;

    }

    public void setIncludeAnnotations(boolean include) {
        myIncludeTaxaAnnotations = include;
    }

    public ExportPlugin setSiteScoreType(SiteScore.SITE_SCORE_TYPE type) {
        mySiteScoreType = type;
        return this;
    }

    public String performFunctionForSimpleTree(SimpleTree input) {

        if (isInteractive()) {
            ReportOptionDialog theDialog = new ReportOptionDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }
            myFileType = theDialog.getTasselFileType();

            theDialog.dispose();

            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String resultFile = Utils.addSuffixIfNeeded(mySaveFile, ".txt");
        if (myFileType == FileLoadPlugin.TasselFileType.Text) {
            BufferedWriter writer = Utils.getBufferedWriter(resultFile);
            try {
                writer.append(input.toString());
            } catch (Exception e) {
                e.printStackTrace();
                throw new IllegalStateException("ExportPlugin: performFunctionForReport: Problem writing file: " + resultFile);
            } finally {
                try {
                    writer.close();
                } catch (Exception e) {
                    // do nothing
                }
            }
        } else {
            PrintWriter writer = null;
            try {
                writer = new PrintWriter(resultFile);
                input.report(writer);
            } catch (Exception e) {
                e.printStackTrace();
                throw new IllegalStateException("ExportPlugin: performFunctionForReport: Problem writing file: " + resultFile);
            } finally {
                try {
                    writer.close();
                } catch (Exception e) {
                    // do nothing
                }
            }
        }
        return resultFile;

    }

    public String performFunctionForTaxaList(TaxaList input) {

        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String filename = "";
        try {
            filename = JSONUtils.exportTaxaListToJSON(input, mySaveFile);
            return new File(filename).getCanonicalPath();
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ExportPlugin: performFunctionForTaxaList: Problem writing file: " + filename);
        }

    }

    public String performFunctionForPositionList(PositionList input) {

        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String filename = "";
        try {
            filename = JSONUtils.exportPositionListToJSON(input, mySaveFile);
            return new File(filename).getCanonicalPath();
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ExportPlugin: performFunctionForPositionList: Problem writing file: " + filename);
        }

    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = ExportPlugin.class.getResource("/net/maizegenetics/analysis/images/Export16.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Export";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Export data to files on your computer.";
    }

    public String getSaveFile() {
        return mySaveFile;
    }

    public ExportPlugin setSaveFile(String saveFile) {
        mySaveFile = saveFile;
        return this;
    }

    public void setSaveFile(File saveFile) {

        if (saveFile == null) {
            mySaveFile = null;
        } else {
            mySaveFile = saveFile.getPath();
        }

    }

    public ExportPlugin setAlignmentFileType(FileLoadPlugin.TasselFileType type) {
        myFileType = type;
        return this;
    }

    private File getFileByChooser() {
        File result = null;
        int returnVal = myFileChooserSave.showSaveDialog(getParentFrame());
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            result = myFileChooserSave.getSelectedFile();
            TasselPrefs.putSaveDir(myFileChooserSave.getCurrentDirectory().getPath());
        } else {
            return null;
        }
        if (result.exists()) {
            int val = JOptionPane.showConfirmDialog(getParentFrame(), "This file already exists: " + result.getName() + "\nDo you want to overwrite it?", "Warning", JOptionPane.YES_NO_OPTION);
            if (val == JOptionPane.NO_OPTION) {
                return null;
            }
        }
        return result;
    }

    class ExportPluginDialog extends JDialog {

        private boolean myIsCancel = true;
        private ButtonGroup myButtonGroup = new ButtonGroup();
        private JRadioButton myHapMapRadioButton = new JRadioButton("Write Hapmap");
        private JRadioButton myByteHDF5RadioButton = new JRadioButton("Write HDF5");
        private JRadioButton myVCFRadioButton = new JRadioButton("Write VCF");
        private JRadioButton myPlinkRadioButton = new JRadioButton("Write Plink");
        private JRadioButton myPhylipRadioButton = new JRadioButton("Write Phylip (Sequential)");
        private JRadioButton myPhylipInterRadioButton = new JRadioButton("Write Phylip (Interleaved)");
        private JRadioButton myTabTableRadioButton = new JRadioButton("Write Tab Delimited");

        private JRadioButton myReferenceProbabilityRadioButton = new JRadioButton("Reference Probability");

        private JCheckBox myKeepDepthCheck = new JCheckBox("Keep Depth (VCF or HDF5)", true);

        private java.util.List<GenotypeTable.GENOTYPE_TABLE_COMPONENT> myComponents;

        public ExportPluginDialog(java.util.List<GenotypeTable.GENOTYPE_TABLE_COMPONENT> components) {
            super((Frame) null, "Export...", true);
            myComponents = components;
            try {
                jbInit();
                pack();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        private void jbInit() throws Exception {

            setTitle("Export...");
            setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
            setUndecorated(false);
            getRootPane().setWindowDecorationStyle(JRootPane.NONE);
            Container contentPane = getContentPane();
            BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
            contentPane.setLayout(layout);
            JPanel main = getMain();
            contentPane.add(main);
            pack();
            setResizable(false);

            myButtonGroup.add(myHapMapRadioButton);
            myButtonGroup.add(myByteHDF5RadioButton);
            myButtonGroup.add(myVCFRadioButton);
            myButtonGroup.add(myPlinkRadioButton);
            myButtonGroup.add(myPhylipRadioButton);
            myButtonGroup.add(myPhylipInterRadioButton);
            myButtonGroup.add(myTabTableRadioButton);

            myButtonGroup.add(myReferenceProbabilityRadioButton);

        }

        private JPanel getMain() {
            JPanel inputs = new JPanel();
            inputs.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
            BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
            inputs.setLayout(layout);
            inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            inputs.add(getLabel());
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            inputs.add(getFileTypePanel());
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            inputs.add(getOptionPanel());
            inputs.add(Box.createRigidArea(new Dimension(1, 5)));
            inputs.add(getButtons());
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            return inputs;
        }

        private JPanel getLabel() {
            JPanel result = new JPanel();
            BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
            result.setLayout(layout);
            result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            JLabel jLabel1 = new JLabel("Choose File Type to Export.");
            jLabel1.setFont(new Font("Dialog", Font.BOLD, 18));
            result.add(jLabel1);
            return result;
        }

        private JPanel getFileTypePanel() {
            JPanel result = new JPanel();
            BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
            result.setLayout(layout);
            result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            result.setBorder(BorderFactory.createEtchedBorder());

            boolean defaultButtonNeedSelected = true;

            result.add(Box.createRigidArea(new Dimension(1, 10)));
            if (myComponents.contains(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype)) {
                result.add(myHapMapRadioButton);
                result.add(myByteHDF5RadioButton);
                result.add(myVCFRadioButton);
                result.add(myPlinkRadioButton);
                result.add(myPhylipRadioButton);
                result.add(myPhylipInterRadioButton);
                result.add(myTabTableRadioButton);

                result.add(Box.createRigidArea(new Dimension(1, 10)));
                myHapMapRadioButton.setSelected(true);
                defaultButtonNeedSelected = false;
            }

            if (!myComponents.isEmpty()) {
                if (myComponents.contains(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability)) {
                    result.add(myReferenceProbabilityRadioButton);
                    if (defaultButtonNeedSelected) {
                        myReferenceProbabilityRadioButton.setSelected(true);
                        defaultButtonNeedSelected = false;
                    }
                }
                result.add(Box.createRigidArea(new Dimension(1, 10)));

            }

            return result;

        }

        private JPanel getOptionPanel() {
            JPanel result = new JPanel();
            BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
            result.setLayout(layout);
            result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            result.setBorder(BorderFactory.createEtchedBorder());
            result.add(myKeepDepthCheck);
            result.add(Box.createRigidArea(new Dimension(1, 10)));

            return result;

        }

        private JPanel getButtons() {

            JButton okButton = new JButton();
            JButton cancelButton = new JButton();

            cancelButton.setText("Cancel");
            cancelButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    cancelButton_actionPerformed(e);
                }
            });

            okButton.setText("OK");
            okButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    okButton_actionPerformed(e);
                }
            });

            JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

            result.add(okButton);

            result.add(cancelButton);

            return result;

        }

        public FileLoadPlugin.TasselFileType getTasselFileType() {
            if (myHapMapRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Hapmap;
            }
            if (myByteHDF5RadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.HDF5;
            }
            if (myVCFRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.VCF;
            }
            if (myPlinkRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Plink;
            }
            if (myPhylipRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Phylip_Seq;
            }
            if (myPhylipInterRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Phylip_Inter;
            }
            if (myTabTableRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Table;
            }
            return null;
        }

        public SiteScore.SITE_SCORE_TYPE getSiteScoreType() {
            if (myReferenceProbabilityRadioButton.isSelected()) {
                return SiteScore.SITE_SCORE_TYPE.ReferenceProbablity;
            }
            return null;
        }

        public boolean keepDepth() {
            return myKeepDepthCheck.isSelected();
        }

        private void okButton_actionPerformed(ActionEvent e) {
            myIsCancel = false;
            setVisible(false);
        }

        private void cancelButton_actionPerformed(ActionEvent e) {
            myIsCancel = true;
            setVisible(false);
        }

        public boolean isCancel() {
            return myIsCancel;
        }
    }
}

class HapmapOptionDialog extends JDialog {

    private boolean exportDiploids = TasselPrefs.getExportPluginExportDiploids();
    private boolean includeTaxaAnnotations = TasselPrefs.getExportPluginIncludeTaxaAnnotations();
    private boolean isCancel = false;
    private final JPanel mainPanel = new JPanel();
    private final JCheckBox diploidCheckBox = new JCheckBox("Export as Diploids");
    private final JCheckBox taxaAnnotationsCheckBox = new JCheckBox("Include Taxa Annotations");
    private final JButton okButton = new JButton("Ok");
    private final JButton cancelButton = new JButton("Cancel");

    public HapmapOptionDialog() {
        super((Frame) null, "Hapmap Options", true);
        initUI();
    }

    private void initUI() {

        setLayout(new BorderLayout());
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
        mainPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        okButton.addActionListener((ActionEvent e) -> {
            okButton_actionPerformed(e);
        });

        cancelButton.addActionListener((ActionEvent e) -> {
            cancelButton_actionPerformed(e);
        });

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);

        diploidCheckBox.setSelected(exportDiploids);
        taxaAnnotationsCheckBox.setSelected(includeTaxaAnnotations);
        mainPanel.add(diploidCheckBox);
        mainPanel.add(taxaAnnotationsCheckBox);

        add(mainPanel, BorderLayout.CENTER);
        add(buttonPanel, BorderLayout.SOUTH);
        pack();
    }

    private void okButton_actionPerformed(ActionEvent e) {
        exportDiploids = diploidCheckBox.isSelected();
        TasselPrefs.putExportPluginExportDiploids(exportDiploids);
        includeTaxaAnnotations = taxaAnnotationsCheckBox.isSelected();
        TasselPrefs.putExportPluginIncludeTaxaAnnotations(includeTaxaAnnotations);
        isCancel = false;
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        isCancel = true;
        setVisible(false);
    }

    public boolean getDiploid() {
        return exportDiploids;
    }

    public boolean includeTaxaAnnotations() {
        return includeTaxaAnnotations;
    }

    public boolean isCancel() {
        return isCancel;
    }
}

class ReportOptionDialog extends JDialog {

    private boolean myIsCancel = true;
    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JRadioButton myReportRadioButton = new JRadioButton("Write As Report");
    private JRadioButton myTextRadioButton = new JRadioButton("Write As Text");

    public ReportOptionDialog() {
        super((Frame) null, "Export Report...", true);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        setTitle("Export Report...");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);

        Container contentPane = getContentPane();

        BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
        contentPane.setLayout(layout);

        JPanel main = getMain();

        contentPane.add(main);

        pack();

        setResizable(false);

        myButtonGroup.add(myReportRadioButton);
        myButtonGroup.add(myTextRadioButton);
        myReportRadioButton.setSelected(true);

    }

    private JPanel getMain() {

        JPanel inputs = new JPanel();
        BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
        inputs.setLayout(layout);
        inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getLabel());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getOptionPanel());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getButtons());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        return inputs;

    }

    private JPanel getLabel() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);

        JLabel jLabel1 = new JLabel("Choose File Type to Export.");
        jLabel1.setFont(new Font("Dialog", Font.BOLD, 18));
        result.add(jLabel1);

        return result;

    }

    private JPanel getOptionPanel() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        result.setBorder(BorderFactory.createEtchedBorder());

        result.add(myReportRadioButton);
        result.add(myTextRadioButton);

        result.add(Box.createRigidArea(new Dimension(1, 20)));

        return result;

    }

    private JPanel getButtons() {

        JButton okButton = new JButton();
        JButton cancelButton = new JButton();

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });

        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });

        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

        result.add(okButton);

        result.add(cancelButton);

        return result;

    }

    public FileLoadPlugin.TasselFileType getTasselFileType() {
        if (myTextRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Text;
        }
        return null;
    }

    private void okButton_actionPerformed(ActionEvent e) {
        myIsCancel = false;
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        myIsCancel = true;
        setVisible(false);
    }

    public boolean isCancel() {
        return myIsCancel;
    }
}

class ExportSquareMatrixDialog extends JDialog {

    private boolean myIsCancel = true;
    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JRadioButton mySquareMatrixButton = new JRadioButton("Write Square Matrix");
    private JRadioButton myRawMultiBlupButton = new JRadioButton("Write Raw MultiBLUP Matrix");
    private JRadioButton myBinMultiBlupButton = new JRadioButton("Write Binary MultiBLUP Matrix");

    public ExportSquareMatrixDialog() {
        super((Frame) null, "Export...", true);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        setTitle("Export...");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);
        Container contentPane = getContentPane();
        BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
        contentPane.setLayout(layout);
        JPanel main = getMain();
        contentPane.add(main);
        pack();
        setResizable(false);

        myButtonGroup.add(mySquareMatrixButton);
        myButtonGroup.add(myRawMultiBlupButton);
        myButtonGroup.add(myBinMultiBlupButton);

        mySquareMatrixButton.setSelected(true);

    }

    private JPanel getMain() {
        JPanel inputs = new JPanel();
        inputs.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
        inputs.setLayout(layout);
        inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        inputs.add(Box.createRigidArea(new Dimension(1, 10)));
        inputs.add(getLabel());
        inputs.add(Box.createRigidArea(new Dimension(1, 10)));
        inputs.add(getFileTypePanel());
        inputs.add(Box.createRigidArea(new Dimension(1, 10)));
        inputs.add(getButtons());
        inputs.add(Box.createRigidArea(new Dimension(1, 10)));
        return inputs;
    }

    private JPanel getLabel() {
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        JLabel jLabel1 = new JLabel("Choose File Type to Export.");
        jLabel1.setFont(new Font("Dialog", Font.BOLD, 18));
        result.add(jLabel1);
        return result;
    }

    private JPanel getFileTypePanel() {
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        result.setBorder(BorderFactory.createEtchedBorder());

        result.add(Box.createRigidArea(new Dimension(1, 10)));
        result.add(mySquareMatrixButton);
        result.add(myRawMultiBlupButton);
        result.add(myBinMultiBlupButton);

        return result;

    }

    private JPanel getButtons() {

        JButton okButton = new JButton();
        JButton cancelButton = new JButton();

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });

        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });

        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

        result.add(okButton);

        result.add(cancelButton);

        return result;

    }

    public FileLoadPlugin.TasselFileType getTasselFileType() {
        if (mySquareMatrixButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.SqrMatrix;
        } else if (myRawMultiBlupButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.SqrMatrixRaw;
        } else if (myBinMultiBlupButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.SqrMatrixBin;
        }
        return FileLoadPlugin.TasselFileType.SqrMatrix;
    }

    private void okButton_actionPerformed(ActionEvent e) {
        myIsCancel = false;
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        myIsCancel = true;
        setVisible(false);
    }

    public boolean isCancel() {
        return myIsCancel;
    }
}
