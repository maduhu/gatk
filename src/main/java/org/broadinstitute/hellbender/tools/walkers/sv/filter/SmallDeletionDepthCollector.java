package org.broadinstitute.hellbender.tools.walkers.sv.filter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.List;

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
}
