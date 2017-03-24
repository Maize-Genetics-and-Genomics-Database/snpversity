/*
 *  PhenotypeUtils
 * 
 *  Created on Oct 27, 2014
 */
package net.maizegenetics.phenotype;

import java.io.BufferedWriter;

import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class PhenotypeUtils {

    private static final Logger myLogger = Logger.getLogger(PhenotypeUtils.class);

    private static final String DELIMITER = "\t";

    private PhenotypeUtils() {
        // utility
    }

    public static void write(Phenotype phenotype, String filename) {

        BufferedWriter writer = null;
        try {
            writer = Utils.getBufferedWriter(filename);

            writer.append("<Phenotype>\n");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if (i != 0) {
                    writer.write(DELIMITER);
                }
                writer.write(phenotype.attributeType(i).name());
            }
            writer.write("\n");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if (i != 0) {
                    writer.write(DELIMITER);
                }
                writer.write(phenotype.attributeName(i));
            }
            writer.write("\n");

            TableReportUtils.saveDelimitedTableReport(phenotype, DELIMITER, writer, false);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("PhenotypeUtils: write: problem saving file: " + filename);
        } finally {
            if (writer != null) {
                try {
                    writer.close();
                } catch (Exception ex) {
                    // do nothing
                }
            }
        }
    }

}
