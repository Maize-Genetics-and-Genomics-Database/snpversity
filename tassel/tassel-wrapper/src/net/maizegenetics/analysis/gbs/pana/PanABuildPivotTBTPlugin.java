package net.maizegenetics.analysis.gbs.pana;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.Arrays;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import static net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups.encodeBySign;
import net.maizegenetics.util.Tassel5HDF5Constants;

/** 
 * Split large TagsByTaxaByteHDF5TagGroup file into small sub TBTs. Designed to submit genetic mapping jobs in cluster
 * 
 * @author Fei Lu
 */
public class PanABuildPivotTBTPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(PanABuildPivotTBTPlugin.class);
    
    String masterTagCountFileS = null;
    String tagCountDirS = null;
    String tbtFileS = null;
    
    IHDF5Writer h5 = null;

    public PanABuildPivotTBTPlugin() {
        super(null, false);
    }

    public PanABuildPivotTBTPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -m  master TagCount file\n"
                + " -d  directory containing tagCount files\n"  
                + " -o  output TBT\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        File[] tcFiles = new File (this.tagCountDirS).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith("cnt");
            }
        });
        Arrays.sort(tcFiles);
        String[] taxaNames = new String[tcFiles.length];
        for (int i = 0; i < taxaNames.length; i++) {
            taxaNames[i] = tcFiles[i].getName().substring(0, tcFiles[i].getName().length()-4);
        }
        TagCounts masterTc = new TagCounts(this.masterTagCountFileS, FilePacking.Byte);
        TagCounts[] tcs = new TagCounts[taxaNames.length];
        for (int i = 0; i < tcs.length; i++) {
            tcs[i] = new TagCounts(tcFiles[i].getAbsolutePath(), FilePacking.Byte);
        }
        this.naiveBuildTBT(masterTc, taxaNames, tcs);
        return null;
    }

    private void naiveBuildTBT (TagCounts masterTC, String[] taxaNames, TagCounts[] tcs) {
        long[][] newTags = new long[masterTC.getTagSizeInLong()][masterTC.getTagCount()];
        for (int i = 0; i < masterTC.getTagCount(); i++) {
            long[] ct = masterTC.getTag(i);
            for (int j = 0; j < masterTC.getTagSizeInLong(); j++) {
                newTags[j][i] = ct[j];
            }
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(this.tbtFileS));
        System.out.println("Creating HDF5 file: " + tbtFileS);
        config.overwrite();
        config.dontUseExtendableDataTypes();
        config.useUTF8CharacterEncoding();
        h5 = config.writer();
        h5.setIntAttribute("/", "tagCount", masterTC.getTagCount());
        h5.setIntAttribute("/", "chunkSize", Tassel5HDF5Constants.BLOCK_SIZE);
        h5.setIntAttribute("/", "tagLengthInLong", masterTC.getTagSizeInLong());
        h5.setIntAttribute("/", "taxaNum", taxaNames.length);
        //create tag matrix
        h5.createLongMatrix("tags", masterTC.getTagSizeInLong(), masterTC.getTagCount(), masterTC.getTagSizeInLong(), masterTC.getTagCount());
        h5.writeLongMatrix("tags", newTags);
        h5.createByteArray("tagLength", masterTC.getTagCount());
        h5.writeByteArray("tagLength", masterTC.getTagLength());
         //create TBT matrix
        h5.createGroup("tbttg");
        int tagChunks = masterTC.getTagCount() >> 16;
        if (masterTC.getTagCount() % Tassel5HDF5Constants.BLOCK_SIZE > 0) {
            tagChunks++;
        }
        System.out.println(Tassel5HDF5Constants.BLOCK_SIZE);
        System.out.printf("tagChunks %d Div %g %n", tagChunks, (double) masterTC.getTagCount() / (double) Tassel5HDF5Constants.BLOCK_SIZE);
        h5.setIntAttribute("tbttg/", "tagCount", masterTC.getTagCount());
        h5.setIntAttribute("tbttg/", "tagChunks", tagChunks);

        for (int tc = 0; tc < tagChunks; tc++) {
            h5.createGroup("tbttg/c" + tc);
        }
        h5.createStringVariableLengthArray("tbttg/taxaNames", taxaNames.length);
        h5.writeStringVariableLengthArray("tbttg/taxaNames", taxaNames);
        
        long[] t;
        byte cnt = 0;
        int index;
        for (int i = 0; i < masterTC.getTagCount(); i++) {
            byte[] td = new byte[taxaNames.length];
            t = masterTC.getTag(i);
            for (int j = 0; j < taxaNames.length; j++) {
                index = tcs[j].getTagIndex(t);
                if (index < 0) cnt = 0;
                else cnt = (byte)tcs[j].getReadCount(index);
                td[j] = cnt;
            }
            int chunk = i >> 16;
            String d = "tbttg/c" + chunk + "/" + i;
            byte[] deftc = encodeBySign(td);
            h5.createByteArray(d, deftc.length);
            h5.writeByteArray(d, deftc);
        }   
    }
    
    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-m", "--master-TC", true);
            engine.add("-d", "--directory-TC", true);
            engine.add("-o", "--output-TBT", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-m")) {
            this.masterTagCountFileS = engine.getString("-m");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-d")) {
            this.tagCountDirS = engine.getString("-d");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        
        if (engine.getBoolean("-o")) {
            this.tbtFileS = engine.getString("-o");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
