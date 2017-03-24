/*
 * FileLoadPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.dna.snp.ReadSequenceAlignmentUtils;
import net.maizegenetics.dna.snp.io.ReadNumericMarkerUtils;
import net.maizegenetics.dna.snp.io.BuilderFromHapMapLIX;
import net.maizegenetics.dna.snp.io.LineIndexBuilder;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.util.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginListener;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.dna.map.TOPMUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.dna.snp.io.FilterJSONUtils;

import org.apache.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.URL;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class FileLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FileLoadPlugin.class);
    private String[] myOpenFiles = null;
    private TasselFileType myFileType = TasselFileType.Unknown;
    private PlinkLoadPlugin myPlinkLoadPlugin = null;
    private ProjectionLoadPlugin myProjectionLoadPlugin = null;
    private ProjectPcsAndRunModelSelectionPlugin myProjectPcsAndRunModelSelectionPlugin = null;
    private final JFileChooser myOpenFileChooser;

    public enum TasselFileType {

        SqrMatrix, Sequence, Unknown, Fasta, Hapmap, HapmapLIX,
        Plink, Phenotype, ProjectionAlignment, ProjectPCsandRunModelSelection, Phylip_Seq, Phylip_Inter, Table,
        Serial, HapmapDiploid, Text, VCF, HDF5, TOPM, HDF5Schema, Filter, NumericGenotype, TaxaList, PositionList, SqrMatrixRaw, SqrMatrixBin
    };
    public static final String FILE_EXT_HAPMAP = ".hmp.txt";
    public static final String FILE_EXT_HAPMAP_GZ = ".hmp.txt.gz";
    public static final String FILE_EXT_HAPMAP_GZ_LIX = FILE_EXT_HAPMAP_GZ + LineIndexBuilder.LINE_INDEX_FILE_EXTENSION;
    public static final String FILE_EXT_PLINK_MAP = ".plk.map";
    public static final String FILE_EXT_PLINK_PED = ".plk.ped";
    public static final String FILE_EXT_SERIAL_GZ = ".serial.gz";
    public static final String FILE_EXT_HDF5 = ".h5";
    public static final String FILE_EXT_VCF = ".vcf";
    public static final String FILE_EXT_TOPM = ".topm";
    public static final String FILE_EXT_TOPM_H5 = ".topm.h5";
    public static final String FILE_EXT_TOPM_BIN = ".topm.bin";
    public static final String FILE_EXT_TOPM_TEXT = ".topm.txt";
    public static final String FILE_EXT_FASTA = ".fasta";

    /**
     * Creates a new instance of FileLoadPlugin
     */
    public FileLoadPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
        if (isInteractive) {
            myOpenFileChooser = new JFileChooser(TasselPrefs.getOpenDir());
            myOpenFileChooser.setMultiSelectionEnabled(true);
        } else {
            myOpenFileChooser = null;
        }
    }

    public DataSet performFunction(DataSet input) {

        myWasCancelled = true;

        try {

            if (isInteractive()) {
                FileLoadPluginDialog theDialog = new FileLoadPluginDialog();
                theDialog.setLocationRelativeTo(getParentFrame());
                theDialog.setVisible(true);
                if (theDialog.isCancel()) {
                    return null;
                }
                myFileType = theDialog.getTasselFileType();

                if (myFileType == TasselFileType.Plink) {
                    if (myPlinkLoadPlugin == null) {
                        myPlinkLoadPlugin = new PlinkLoadPlugin(getParentFrame(), isInteractive());
                        for (PluginListener current : getListeners()) {
                            myPlinkLoadPlugin.addListener(current);
                        }
                    }
                    return myPlinkLoadPlugin.performFunction(null);
                }

                if (myFileType == TasselFileType.ProjectionAlignment) {
                    if (myProjectionLoadPlugin == null) {
                        myProjectionLoadPlugin = new ProjectionLoadPlugin(getParentFrame(), isInteractive());
                        for (PluginListener current : getListeners()) {
                            myProjectionLoadPlugin.addListener(current);
                        }
                    }
                    return myProjectionLoadPlugin.performFunction(input);
                }

                if (myFileType == TasselFileType.ProjectPCsandRunModelSelection) {
                    if (myProjectPcsAndRunModelSelectionPlugin == null) {
                        myProjectPcsAndRunModelSelectionPlugin = new ProjectPcsAndRunModelSelectionPlugin(getParentFrame(), isInteractive());
                        for (PluginListener current : getListeners()) {
                            myProjectPcsAndRunModelSelectionPlugin.addListener(current);
                        }
                    }
                    return myProjectPcsAndRunModelSelectionPlugin.performFunction(input);
                }
                setOpenFiles(getOpenFilesByChooser());
                theDialog.dispose();
            }

            if ((myOpenFiles == null) || (myOpenFiles.length == 0)) {
                return null;
            }

            List result = new ArrayList();
            ArrayList<String> alreadyLoaded = new ArrayList();
            for (int i = 0; i < myOpenFiles.length; i++) {

                if (alreadyLoaded.contains(myOpenFiles[i])) {
                    continue;
                }

                LocalDateTime time = LocalDateTime.now();
                String timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
                myLogger.info("Start Loading File: " + myOpenFiles[i] + " time: " + timeStr);

                DataSet tds = null;
                try {

                    if (myFileType == TasselFileType.Unknown) {
                        if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP_GZ)) {
                            String theIndex = myOpenFiles[i].replaceFirst(FILE_EXT_HAPMAP_GZ, FILE_EXT_HAPMAP_GZ_LIX);
                            if (new File(theIndex).isFile()) {
                                myLogger.info("guessAtUnknowns: type: " + TasselFileType.HapmapLIX);
                                alreadyLoaded.add(myOpenFiles[i]);
                                alreadyLoaded.add(theIndex);
                                GenotypeTable hapmap = BuilderFromHapMapLIX.build(myOpenFiles[i], theIndex);
                                tds = new DataSet(new Datum(Utils.getFilename(myOpenFiles[i], FileLoadPlugin.FILE_EXT_HAPMAP_GZ), hapmap, null), this);
                            } else {
                                myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                                alreadyLoaded.add(myOpenFiles[i]);
                                tds = processDatum(myOpenFiles[i], TasselFileType.Hapmap);
                            }
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP_GZ_LIX)) {
                            String theHapmap = myOpenFiles[i].replaceFirst(FILE_EXT_HAPMAP_GZ_LIX, FILE_EXT_HAPMAP_GZ);
                            if (new File(theHapmap).isFile()) {
                                myLogger.info("guessAtUnknowns: type: " + TasselFileType.HapmapLIX);
                                alreadyLoaded.add(myOpenFiles[i]);
                                alreadyLoaded.add(theHapmap);
                                GenotypeTable hapmap = BuilderFromHapMapLIX.build(theHapmap, myOpenFiles[i]);
                                tds = new DataSet(new Datum(Utils.getFilename(theHapmap, FileLoadPlugin.FILE_EXT_HAPMAP_GZ), hapmap, null), this);
                            } else {
                                throw new IllegalStateException("FileLoadPlugin: Can't find file matching: " + myOpenFiles[i]);
                            }
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Hapmap);
                        } else if ((myOpenFiles[i].endsWith(FILE_EXT_TOPM_H5)) || (myOpenFiles[i].endsWith(FILE_EXT_TOPM))
                                || (myOpenFiles[i].endsWith(FILE_EXT_TOPM_BIN)) || (myOpenFiles[i].endsWith(FILE_EXT_TOPM_TEXT))) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.TOPM);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.TOPM);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_PED)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                            String theMapFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_PED, FILE_EXT_PLINK_MAP);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(theMapFile);
                            GenotypeTable plink = ImportUtils.readFromPLink(myOpenFiles[i], theMapFile, this);
                            tds = new DataSet(new Datum(Utils.getFilename(myOpenFiles[i], FileLoadPlugin.FILE_EXT_PLINK_PED), plink, null), this);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_MAP)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                            String thePedFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_MAP, FILE_EXT_PLINK_PED);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(thePedFile);
                            GenotypeTable plink = ImportUtils.readFromPLink(thePedFile, myOpenFiles[i], this);
                            tds = new DataSet(new Datum(Utils.getFilename(thePedFile, FileLoadPlugin.FILE_EXT_PLINK_PED), plink, null), this);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_SERIAL_GZ)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Serial);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Serial);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_HDF5)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.HDF5);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.HDF5);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_VCF) || myOpenFiles[i].endsWith(FILE_EXT_VCF + ".gz")) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.VCF);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.VCF);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_FASTA) || myOpenFiles[i].endsWith(FILE_EXT_FASTA + ".gz")) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Fasta);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Fasta);
                        } else {
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = guessAtUnknowns(myOpenFiles[i]);
                        }
                    } else {
                        alreadyLoaded.add(myOpenFiles[i]);
                        tds = processDatum(myOpenFiles[i], myFileType);
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    StringBuilder builder = new StringBuilder();
                    builder.append("Error loading: ");
                    builder.append(myOpenFiles[i]);
                    builder.append("\n");
                    builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
                    String str = builder.toString();
                    if (isInteractive()) {
                        DialogUtils.showError(str, getParentFrame());
                    } else {
                        myLogger.error(str);
                    }
                }

                time = LocalDateTime.now();
                timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
                if (tds != null) {
                    myLogger.info("Finished Loading File: " + myOpenFiles[i] + " time: " + timeStr);
                    GenotypeSummaryPlugin.printSimpleSummary(tds);
                    myWasCancelled = false;
                    result.add(tds);
                    fireDataSetReturned(new PluginEvent(tds, FileLoadPlugin.class));
                } else {
                    myLogger.info("Nothing Loaded for File: " + myOpenFiles[i] + " time: " + timeStr);
                }

            }

            return DataSet.getDataSet(result, this);

        } catch (Exception e) {
            showError(e, null);
            return null;
        } finally {
            fireProgress(100);
        }

    }

    public DataSet guessAtUnknowns(String filename) {

        TasselFileType guess = TasselFileType.Sequence;
        DataSet tds = null;

        try (BufferedReader br = Utils.getBufferedReader(filename)) {

            String line1 = br.readLine();
            while (line1 != null) {
                line1 = line1.trim();
                if (!line1.isEmpty()) {
                    break;
                }
                line1 = br.readLine();
            }
            if (line1 == null) {
                throw new IllegalArgumentException("FileLoadPlugin: guessAtUnknowns: File is empty: " + filename);
            }
            String[] sval1 = line1.split("\\s");
            String line2 = br.readLine().trim();
            String[] sval2 = line2.split("\\s");
            if (line1.startsWith("{")) {
                String temp;
                if (sval1.length > 1) {
                    temp = sval1[1];
                } else {
                    temp = line2;
                }
                if (temp.startsWith("\"TaxaList\"")) {
                    guess = TasselFileType.TaxaList;
                } else if (temp.startsWith("\"PositionList\"")) {
                    guess = TasselFileType.PositionList;
                } else if (temp.startsWith("\"Filter\"")) {
                    guess = TasselFileType.Filter;
                }
            } else if (line1.startsWith("<") || line1.startsWith("#")) {
                boolean isTrait = false;
                boolean isMarker = false;
                boolean isNumeric = false;
                Pattern tagPattern = Pattern.compile("[<>\\s]+");
                String[] info1 = tagPattern.split(line1);
                String[] info2 = tagPattern.split(line2);
                if (info1.length > 1) {
                    if (info1[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info1[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info1[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    } else if (info1[1].toUpperCase().startsWith("PHENO")) {
                        isTrait = true;
                    }
                }
                if (info2.length > 1) {
                    if (info2[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info2[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info2[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    }
                } else {
                    guess = null;
                    String inline = br.readLine();
                    while (guess == null && inline != null && (inline.startsWith("#") || inline.startsWith("<"))) {
                        if (inline.startsWith("<")) {
                            String[] info = tagPattern.split(inline);
                            if (info[1].toUpperCase().startsWith("MARKER")) {
                                isMarker = true;
                            } else if (info[1].toUpperCase().startsWith("TRAIT")) {
                                isTrait = true;
                            } else if (info[1].toUpperCase().startsWith("NUMER")) {
                                isNumeric = true;
                            }
                        }
                    }
                }
                if (isTrait) {
                    guess = TasselFileType.Phenotype;
                } else if (isMarker && isNumeric) {
                    guess = TasselFileType.NumericGenotype;
                } else {
                    myLogger.warn("Line1: " + line1);
                    myLogger.warn("Line2: " + line2);
                    throw new IOException("Improperly formatted header. Data will not be imported.");
                }
            } else if ((line1.startsWith(">")) || (line1.startsWith(";"))) {
                guess = TasselFileType.Fasta;
            } else if (sval1.length == 1) {
                guess = TasselFileType.SqrMatrix;
            } else if ((line1.startsWith("#Nexus")) || (line1.startsWith("#NEXUS")) || (line1.startsWith("CLUSTAL"))
                    || ((sval1.length == 2) && (sval2.length == 2))) {
                guess = TasselFileType.Sequence;
            }

            myLogger.info("guessAtUnknowns: type: " + guess);
            tds = processDatum(filename, guess);

        } catch (Exception e) {
            showError(e, filename);
        }

        return tds;

    }

    private DataSet processDatum(String inFile, TasselFileType theFT) {
        Object result = null;
        String suffix = null;
        try {
            switch (theFT) {
                case Hapmap: {
                    suffix = FILE_EXT_HAPMAP;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_HAPMAP_GZ;
                    }
                    result = ImportUtils.readFromHapmap(inFile, this);
                    break;
                }
                case HDF5: {
                    IHDF5Reader reader = HDF5Factory.openForReading(inFile);
                    boolean t4HDF5 = HDF5Utils.isTASSEL4HDF5Format(HDF5Factory.openForReading(inFile));
                    reader.close();
                    if (t4HDF5) {
                        String newInfile = inFile.replace(".h5", ".t5.h5");
                        if (new File(newInfile).exists()) {
                            String message = "This file is TASSEL 4 HDF5 format. It looks like it has already been converted to TASSEL 5. Using file: " + newInfile;
                            if (isInteractive()) {
                                DialogUtils.showWarning(message, getParentFrame());
                            } else {
                                myLogger.warn(message);
                            }
                        } else {
                            String message = "This file is TASSEL 4 HDF5 format. It will be converted to TASSEL 5 "
                                    + "HDF5 format with name: " + newInfile + ".  This may take a few minutes.";
                            if (isInteractive()) {
                                DialogUtils.showWarning(message, getParentFrame());
                            } else {
                                myLogger.warn(message);
                            }
                            MigrateHDF5FromT4T5.copyGenotypes(inFile, newInfile);
                        }

                        inFile = newInfile;
                    }
                    suffix = FILE_EXT_HDF5;
                    result = ImportUtils.readGuessFormat(inFile);
                    break;
                }
                case HDF5Schema: {
                    suffix = "";
                    result = new HDF5TableReport(inFile);
                    break;
                }
                case VCF: {
                    suffix = FILE_EXT_VCF;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_VCF + ".gz";
                    }
                    result = ImportUtils.readFromVCF(inFile, this);
                    break;
                }
                case Sequence: {
                    result = ReadSequenceAlignmentUtils.readBasicAlignments(inFile, 40);
                    break;
                }
                case Fasta: {
                    result = ImportUtils.readFasta(inFile);
                    break;
                }
                case SqrMatrix: {
                    result = ReadDistanceMatrix.readDistanceMatrix(inFile);
                    break;
                }
                case Phenotype: {
                    List<Phenotype> phenotypes = new PhenotypeBuilder().fromFile(inFile).build();
                    if (phenotypes.size() != 1) {
                        throw new IllegalStateException("FileLoadPlugin: processDatum: problem loading phenotype file: " + inFile);
                    }
                    result = phenotypes.get(0);
                    break;
                }
                case NumericGenotype: {
                    result = ReadNumericMarkerUtils.readNumericMarkerFile(inFile);
                    break;
                }
                case TaxaList: {
                    result = JSONUtils.importTaxaListFromJSON(inFile);
                    break;
                }
                case PositionList: {
                    result = JSONUtils.importPositionListFromJSON(inFile);
                    break;
                }
                case Table: {
                    result = TableReportUtils.readDelimitedTableReport(inFile, "\t");
                    break;
                }
                case TOPM: {
                    result = TOPMUtils.readTOPM(inFile);
                    break;
                }
                case Filter: {
                    result = FilterJSONUtils.importJSONToFilter(inFile);
                    break;
                }
            }
        } catch (Exception e) {
            showError(e, inFile);
        }
        if (result != null) {
            String name = Utils.getFilename(inFile, suffix);

            Datum td = new Datum(name, result, null);
            //todo need to add logic of directories.
            DataSet tds = new DataSet(td, this);
            return tds;
        }
        return null;
    }

    private void showError(Exception e, String filename) {

        myLogger.error(e.getMessage(), e);
        StringBuilder builder = new StringBuilder();
        if ((filename != null) && (filename.length() != 0)) {
            builder.append("Error loading: ");
            builder.append(filename);
            builder.append("\n");
        }
        builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
        String str = builder.toString();
        if (isInteractive()) {
            DialogUtils.showError(str, getParentFrame());
        } else {
            myLogger.error(str);
        }

    }

    /**
     * Provides a open filer that remember the last location something was
     * opened from
     */
    private File[] getOpenFilesByChooser() {
        File[] lopenFiles = null;
        myOpenFileChooser.setVisible(true);
        int returnVal = myOpenFileChooser.showOpenDialog(getParentFrame());
        if (returnVal == JFileChooser.OPEN_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            lopenFiles = myOpenFileChooser.getSelectedFiles();
            TasselPrefs.putOpenDir(myOpenFileChooser.getCurrentDirectory().getPath());
        }
        return lopenFiles;
    }

    public String[] getOpenFiles() {
        return myOpenFiles;
    }

    public void setOpenFiles(File[] openFiles) {

        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
            return;
        }

        myOpenFiles = new String[openFiles.length];
        for (int i = 0; i < openFiles.length; i++) {
            myOpenFiles[i] = openFiles[i].getPath();
        }

    }

    public void setOpenFiles(String[] openFiles) {
        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
        } else {
            myOpenFiles = openFiles;
        }
    }

    public TasselFileType getTheFileType() {
        return myFileType;
    }

    public void setTheFileType(TasselFileType theFileType) {
        myFileType = theFileType;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FileLoadPlugin.class.getResource("/net/maizegenetics/analysis/images/LoadFile.gif");
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
        return "Load";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load data from files on your computer.";
    }
}

/**
 * <p>Title: TASSEL</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p>Company: USDA-ARS</p>
 *
 * @author Edward Buckler
 * @version 1.0
 */
class FileLoadPluginDialog extends JDialog {

    boolean isCancel = true;
    ButtonGroup conversionButtonGroup = new ButtonGroup();
    JRadioButton hapMapRadioButton = new JRadioButton("Load Hapmap");
    JRadioButton hdf5RadioButton = new JRadioButton("Load HDF5");
    JRadioButton hdf5SchemaRadioButton = new JRadioButton("Load HDF5 Schema");
    JRadioButton vcfRadioButton = new JRadioButton("Load VCF");
    JRadioButton plinkRadioButton = new JRadioButton("Load Plink");
    JRadioButton sequenceAlignRadioButton = new JRadioButton("Load Phylip");
    JRadioButton fastaRadioButton = new JRadioButton("Load FASTA File");
    JRadioButton numericalRadioButton = new JRadioButton("Load Numerical (trait, covariates, or factors)");
    JRadioButton loadMatrixRadioButton = new JRadioButton("Load Square Numerical Matrix (i.e. kinship)");
    JRadioButton guessRadioButton = new JRadioButton("Make Best Guess");
    JRadioButton projectionAlignmentRadioButton = new JRadioButton("Load Projection Alignment");
    JRadioButton projectPCsandRunModelSelectionRadioButton = new JRadioButton("Load Files for Projecting PCs onto NAM");
    JRadioButton tableReportRadioButton = new JRadioButton("Load a Table Report");
    JRadioButton topmRadioButton = new JRadioButton("Load a TOPM (Tags on Physical Map)");

    public FileLoadPluginDialog() {
        super((Frame) null, "File Loader", true);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        setTitle("File Loader");
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

        conversionButtonGroup.add(projectionAlignmentRadioButton);
        conversionButtonGroup.add(projectPCsandRunModelSelectionRadioButton);
        conversionButtonGroup.add(hapMapRadioButton);
        conversionButtonGroup.add(hdf5RadioButton);
        conversionButtonGroup.add(hdf5SchemaRadioButton);
        conversionButtonGroup.add(vcfRadioButton);
        conversionButtonGroup.add(plinkRadioButton);
        conversionButtonGroup.add(sequenceAlignRadioButton);
        conversionButtonGroup.add(fastaRadioButton);
        conversionButtonGroup.add(loadMatrixRadioButton);
        conversionButtonGroup.add(numericalRadioButton);
        conversionButtonGroup.add(tableReportRadioButton);
        conversionButtonGroup.add(topmRadioButton);
        conversionButtonGroup.add(guessRadioButton);
        guessRadioButton.setSelected(true);

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

        JLabel jLabel1 = new JLabel("Choose File Type to Load.");
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

        result.add(hapMapRadioButton);
        result.add(hdf5RadioButton);
        result.add(hdf5SchemaRadioButton);
        result.add(vcfRadioButton);
        result.add(plinkRadioButton);
        result.add(projectionAlignmentRadioButton);
        //result.add(projectPCsandRunModelSelectionRadioButton);
        result.add(sequenceAlignRadioButton);
        result.add(fastaRadioButton);
        result.add(numericalRadioButton);
        result.add(loadMatrixRadioButton);
        result.add(tableReportRadioButton);
        result.add(topmRadioButton);
        result.add(guessRadioButton);

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
        if (hapMapRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Hapmap;
        }
        if (hdf5RadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.HDF5;
        }
        if (hdf5SchemaRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.HDF5Schema;
        }
        if (vcfRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.VCF;
        }
        if (plinkRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Plink;
        }
        if (projectionAlignmentRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.ProjectionAlignment;
        }
        if (projectPCsandRunModelSelectionRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.ProjectPCsandRunModelSelection;
        }
        if (sequenceAlignRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Sequence;
        }
        if (fastaRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Fasta;
        }
        if (loadMatrixRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.SqrMatrix;
        }
        if (numericalRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Unknown;
        }
        if (tableReportRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Table;
        }
        if (topmRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.TOPM;
        }
        return FileLoadPlugin.TasselFileType.Unknown;
    }

    public void okButton_actionPerformed(ActionEvent e) {
        isCancel = false;
        setVisible(false);
    }

    public void cancelButton_actionPerformed(ActionEvent e) {
        isCancel = true;
        setVisible(false);
    }

    public boolean isCancel() {
        return isCancel;
    }
}
