//
// ReadDistanceMatrix.java
//
package net.maizegenetics.taxa.distance;

import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;

import java.io.IOException;
import java.io.BufferedReader;

import java.util.regex.Pattern;

import org.apache.log4j.Logger;

/**
 * @author Terry Casstevens
 */
public class ReadDistanceMatrix {

    private static final Logger myLogger = Logger.getLogger(ReadDistanceMatrix.class);

    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s+");

    private ReadDistanceMatrix() {
    }

    public static DistanceMatrix readDistanceMatrix(String filename) {

        try (BufferedReader reader = Utils.getBufferedReader(filename)) {

            int numSeqs = Integer.parseInt(reader.readLine().trim());

            double[][] distance = new double[numSeqs][numSeqs];
            TaxaListBuilder taxa = new TaxaListBuilder();
            String current = reader.readLine();
            int index = 0;
            while (current != null) {

                if (index >= numSeqs) {
                    throw new IllegalArgumentException("ReadDistancMatrix: There are too many lines in this file.  Expected: " + (numSeqs + 1) + " counting the first line (number of taxa)");
                }

                String[] tokens = WHITESPACE_PATTERN.split(current);
                if (tokens.length != numSeqs + 1) {
                    throw new IllegalStateException("ReadDistanceMatrix: Incorrect number of values on line number: " + (index + 2) + " expected: " + (numSeqs + 1) + " counting taxon name. actual: " + tokens.length);
                }
                taxa.add(new Taxon(tokens[0]));

                for (int i = 0; i < numSeqs; i++) {
                    try {
                        distance[index][i] = Double.parseDouble(tokens[i + 1]);
                    } catch (NumberFormatException nfex) {
                        myLogger.debug(nfex.getMessage(), nfex);
                        throw new IllegalArgumentException("ReadDistanceMatrix: Incorrectly formatted number: " + tokens[i + 1] + " on line number: " + (index + 2));
                    }
                }

                current = reader.readLine();
                index++;
            }

            if (index != numSeqs) {
                throw new IllegalArgumentException("ReadDistanceMatrix: There are too few lines in this file.  Expected: " + (numSeqs + 1) + " counting the first line (number of taxa)");
            }

            return new DistanceMatrix(distance, taxa.build());

        } catch (IOException ioex) {
            myLogger.debug(ioex.getMessage(), ioex);
            throw new IllegalStateException("ReadDistanceMatrix: Problem reading file: " + filename);
        }

    }
}
