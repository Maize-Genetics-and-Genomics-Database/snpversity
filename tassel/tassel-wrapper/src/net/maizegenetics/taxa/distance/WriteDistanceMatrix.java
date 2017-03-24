package net.maizegenetics.taxa.distance;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 * @author Terry Casstevens
 */
public class WriteDistanceMatrix {

    private static final Logger myLogger = Logger.getLogger(WriteDistanceMatrix.class);

    private WriteDistanceMatrix() {
    }

    public static void saveDelimitedDistanceMatrix(DistanceMatrix matrix, String saveFile) {

        if ((saveFile == null) || (saveFile.isEmpty())) {
            throw new IllegalArgumentException("WriteDistanceMatrix: saveDelimitedDistanceMatrix: No file specified.");
        }

        try (BufferedWriter bw = Utils.getBufferedWriter(saveFile)) {

            bw.write(String.valueOf(matrix.getRowCount()));
            bw.write("\n");

            for (long r = 0, n = matrix.getRowCount(); r < n; r++) {
                Object[] theRow = matrix.getRow(r);
                for (int i = 0; i < theRow.length; i++) {
                    if (i != 0) {
                        bw.write("\t");
                    }
                    bw.write(theRow[i].toString());
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("WriteDistanceMatrix: saveDelimitedDistanceMatrix: problem writing file: " + e.getMessage());
        }

    }

    public static void saveRawMultiBlupMatrix(DistanceMatrix matrix, File taxaFile, File matrixFile) {

        if (matrixFile == null || taxaFile == null) {
            return;
        }
        FileWriter taxafw = null;
        BufferedWriter taxabw = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {
            //Write the taxa/ids.  Need 2 columns for MultiBlup.
            taxafw = new FileWriter(taxaFile);
            taxabw = new BufferedWriter(taxafw);

            //Write the matrix
            fw = new FileWriter(matrixFile);
            bw = new BufferedWriter(fw);

            for (long r = 0, n = matrix.getRowCount(); r < n; r++) {
                Object[] theRow = matrix.getRow(r);
                taxabw.write(theRow[0].toString() + " " + theRow[0].toString() + "\n");

                for (int i = 1; i < theRow.length; i++) {
                    if (i != 1) {
                        bw.write("\t");
                    }
                    bw.write(theRow[i].toString());
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            System.out.println("WriteDistanceMatrix: saveDelimitedDistanceMatrix: problem writing file: " + e.getMessage());
        } finally {
            try {
                taxabw.flush();
                bw.flush();
                taxafw.close();
                taxabw.close();
                bw.close();
                fw.close();
            } catch (Exception e) {
                System.out.println(e);
                // do nothing
            }
        }

    }

    public static void saveBinMultiBlupMatrix(DistanceMatrix matrix, File taxaFile, File matrixFile, File countFile) {

        if (matrixFile == null || taxaFile == null || countFile == null) {
            return;
        }
        FileWriter taxafw = null;
        BufferedWriter taxabw = null;
        FileOutputStream countfw = null;
        BufferedOutputStream countbw = null;
        FileOutputStream fw = null;
        BufferedOutputStream bw = null;
        try {
            //Write the taxa/ids.  Need 2 columns for MultiBlup.
            taxafw = new FileWriter(taxaFile);
            taxabw = new BufferedWriter(taxafw);

            countfw = new FileOutputStream(countFile);
            countbw = new BufferedOutputStream(countfw);

            //Write the matrix
            fw = new FileOutputStream(matrixFile);
            bw = new BufferedOutputStream(fw);

            for (long r = 0, n = matrix.getRowCount(); r < n; r++) {
                Object[] theRow = matrix.getRow(r);

                taxabw.write(theRow[0].toString() + " " + theRow[0].toString() + "\n");
                for (int i = 1; i < r + 2; i++) {

                    //ByteBuffer kinsBuffer = ByteBuffer.allocate(4).order(ByteOrder.nativeOrder());
                    ByteBuffer kinsBuffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);

                    kinsBuffer.putFloat(Float.parseFloat(theRow[i].toString()));
                    bw.write(kinsBuffer.array());

                    //ByteBuffer countsBuffer = ByteBuffer.allocate(4).order(ByteOrder.nativeOrder());
                    ByteBuffer countsBuffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);

                    countsBuffer.putFloat(1f);
                    countbw.write(countsBuffer.array());
                }
            }

        } catch (Exception e) {
            System.out.println("WriteDistanceMatrix: saveDelimitedDistanceMatrix: problem writing file: " + e.getMessage());
        } finally {
            try {
                taxabw.flush();
                countbw.flush();
                bw.flush();
                taxafw.close();
                taxabw.close();
                countfw.close();
                countbw.close();
                bw.close();
                fw.close();
            } catch (Exception e) {
                System.out.println(e);
                // do nothing
            }
        }

    }
}
