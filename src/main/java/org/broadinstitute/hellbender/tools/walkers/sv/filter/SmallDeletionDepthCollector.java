package org.broadinstitute.hellbender.tools.walkers.sv.filter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.bed.SimpleBEDFeature;
import org.apache.logging.log4j.Logger;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;
import scala.Tuple3;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;

public class SmallDeletionDepthCollector {

    // TODO: 7/26/18 should secondary (256 sam flag) reads be filtered as well?
    public static List<ReadFilter> customReadFilter() {
        return Arrays.asList(ReadFilterLibrary.MAPPED,
                ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
    }

    public static void writeBedFile(final String output, final SAMSequenceDictionary sequenceDictionary,
                                    final List<AnnotatedIntervalTrio> annotatedIntervalTrios) {

        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(output)))) {
            for (final AnnotatedIntervalTrio annotatedIntervalTrio : annotatedIntervalTrios) {
                writer.write(annotatedIntervalTrio.toBed_5_String(sequenceDictionary));
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write " + output, ioe);
        }

        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(output.replace(".bed", ".summary"))))) {
            for (final AnnotatedIntervalTrio annotatedIntervalTrio : annotatedIntervalTrios) {
                writer.write(annotatedIntervalTrio.toSummaryString(sequenceDictionary));
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write summary for " + output, ioe);
        }
    }
    
    private static SVInterval makeCStyleInterval(final Locatable locatable, final SAMSequenceDictionary samSequenceDictionary) {
        return new SVInterval(samSequenceDictionary.getSequenceIndex(locatable.getContig()), locatable.getStart() - 1, locatable.getEnd());
    }

    // for walker ======================================================================================================

    static SimpleInterval getEqualSizeFlankedInterval(final SAMSequenceDictionary sequenceDictionary, final SimpleBEDFeature feature) {
        int lengthOnReference = feature.getLengthOnReference();

        // TODO: 7/26/18 this based on a quick scan of htsjdk code that BEDFeature, being in tribble, returns Locatable-like coordinates, i.e. [1, ...]
        final String contig = feature.getContig();
        int start = Math.max(0, feature.getStart() - lengthOnReference);
        int end = Math.min(feature.getEnd() + lengthOnReference, sequenceDictionary.getSequence(contig).getSequenceLength());

        return new SimpleInterval(contig, start, end);
    }

    static AnnotatedIntervalTrio forWalkerApply(final SimpleBEDFeature feature, final ReadsContext readsContext,
                                                final ReferenceContext referenceContext, final FeatureContext featureContext,
                                                final SAMSequenceDictionary sequenceDictionary) {

        final String name = feature.getName();
        final SVInterval featureInterval = SmallDeletionDepthCollector.makeCStyleInterval(feature, sequenceDictionary);
        final DepthAnnotatedCStyleInterval middle = new DepthAnnotatedCStyleInterval(new NamedCStyleInterval(featureInterval, name));
        final DepthAnnotatedCStyleInterval leftFlank = new DepthAnnotatedCStyleInterval(
                new NamedCStyleInterval(middle.getContig(),
                        Math.max(0, middle.getStart() - middle.getLength()),
                        middle.getStart(), name));
        int chrLength = sequenceDictionary.getSequence(middle.getContig()).getSequenceLength();
        final DepthAnnotatedCStyleInterval rightFlank = new DepthAnnotatedCStyleInterval(
                new NamedCStyleInterval(middle.getContig(),
                        middle.getEnd(),
                        Math.min(middle.getEnd() + middle.getLength(), chrLength),
                        name));
        final AnnotatedIntervalTrio annotatedIntervalTrio = new AnnotatedIntervalTrio(leftFlank, middle, rightFlank);
        for (final GATKRead gatkRead : readsContext) {
            final SVInterval readInterval = SmallDeletionDepthCollector.makeCStyleInterval(gatkRead, sequenceDictionary); // this uses the read's clipped start
            leftFlank.collectFromRead(gatkRead, readInterval);
        }
        return annotatedIntervalTrio;
    }

    // for spark =======================================================================================================

    /**
     * Given interval of interest, generate left and right flanks of the same size of the interval.
     */
    public static SVIntervalTree<Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval>>
    makeFlankingTree(final String inputBedFile, final SAMSequenceDictionary sequenceDictionary, final Logger logger) {

        final SVIntervalTree<Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval>> flankedTree = new SVIntervalTree<>();

        final BufferedLineReader reader = new BufferedLineReader(BucketUtils.openFile(inputBedFile));
        String line;
        int linenumber = 0;
        while ( (line = reader.readLine()) != null) {
            ++linenumber;
            final String[] bedLine = line.split("\t");
            if (bedLine.length < 3)
                throw new UserException("Line " + linenumber + " ill-formed: " + line);
            final SVInterval svInterval = new SVInterval(sequenceDictionary.getSequenceIndex(bedLine[0]),
                    Integer.valueOf(bedLine[1]), Integer.valueOf(bedLine[2]));
            final StringBuilder builder = new StringBuilder();
            for (int i = 3; i < bedLine.length; ++i) {
                builder.append(bedLine[i]);
            }
            String name = builder.toString();
            final NamedCStyleInterval middle = new NamedCStyleInterval(svInterval, name);
            int chrLength = sequenceDictionary.getSequence(middle.getContig()).getSequenceLength();
            final NamedCStyleInterval leftFlank = new NamedCStyleInterval(middle.getContig(),
                    Math.max(0, middle.getStart() - middle.getLength()),
                    middle.getStart(),
                    name);
            final NamedCStyleInterval rightFlank = new NamedCStyleInterval(middle.getContig(),
                    middle.getEnd(),
                    Math.min(middle.getEnd() + middle.getLength(), chrLength),
                    name);
            flankedTree.put(new SVInterval(middle.getContig(), leftFlank.getStart(), rightFlank.getEnd()),
                    new Tuple3<>(leftFlank, middle, rightFlank));

        }
        logger.info("processing " + linenumber + " intervals");
        return flankedTree;
    }

    // sometimes the flanked intervals intersect so a read could generate multiple tree nodes, i.e. there's overlap between the padded intervals
    public static Iterator<Tuple2<SVInterval, AnnotatedIntervalTrio>>
    collectSummaries(final GATKRead gatkRead, final Broadcast<SAMSequenceDictionary> sequenceDictionaryBroadcast,
                     final Broadcast<SVIntervalTree<Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval>>> treeBroadcast) {

        final SVInterval readInterval = makeCStyleInterval(gatkRead, sequenceDictionaryBroadcast.getValue()); // this uses the read's clipped start

        final Iterator<SVIntervalTree.Entry<Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval>>> overlappersIterator =
                treeBroadcast.getValue().overlappers(readInterval);
        if (overlappersIterator.hasNext()) {
            final List<Tuple2<SVInterval, AnnotatedIntervalTrio>> result = new ArrayList<>(3); // 3 is a guess
            while (overlappersIterator.hasNext()) {
                final SVIntervalTree.Entry<Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval>> entry = overlappersIterator.next();
                result.add(new Tuple2<>(entry.getInterval(), new AnnotatedIntervalTrio(gatkRead, readInterval, entry.getValue())));
            }
            return result.iterator();
        } else
            return Collections.emptyIterator();
    }

}
