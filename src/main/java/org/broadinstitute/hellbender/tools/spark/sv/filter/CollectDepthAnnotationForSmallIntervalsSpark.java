package org.broadinstitute.hellbender.tools.spark.sv.filter;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.walkers.sv.filter.AnnotatedIntervalTrio;
import org.broadinstitute.hellbender.tools.walkers.sv.filter.NamedCStyleInterval;
import org.broadinstitute.hellbender.tools.walkers.sv.filter.SmallDeletionDepthCollector;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "WARNING: THIS DOES NOT ALWAYS SUCCEED DUE TO A (POSSIBLE) BUG IN HADOOP-BAM",
        summary =
                "This tool packages the algorithms described in FindBreakpointEvidenceSpark and" +
                        " DiscoverVariantsFromContigAlignmentsSAMSpark as an integrated workflow.  Please consult the" +
                        " descriptions of those tools for more details about the algorithms employed.  In brief, input reads are examined" +
                        " for evidence of structural variation in a genomic region, regions so identified are locally assembled, and" +
                        " the local assemblies are called for structural variation.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class CollectDepthAnnotationForSmallIntervalsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "path to write BED file for the intervals and flanking regions",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String bedFile;

    @Argument(doc = "path to input BED file holding intervals of interest (trailing columns will be kept)", shortName = "bed", fullName = "inputBED")
    private String inputBED;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultReadFilters = new ArrayList<>( super.getDefaultReadFilters() );
        defaultReadFilters.addAll( SmallDeletionDepthCollector.customReadFilter() );
        return defaultReadFilters;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final SAMSequenceDictionary sequenceDictionary = getHeaderForReads().getSequenceDictionary();

        // first get flanks
        final SVIntervalTree<Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval>> intervalTree =
                SmallDeletionDepthCollector.makeFlankingTree(inputBED, sequenceDictionary, logger);

        final List<SimpleInterval> paddedIntervals = Utils.stream(intervalTree.iterator())
                .map(SVIntervalTree.Entry::getInterval)
                .map(svInterval -> new SVInterval(svInterval.getContig(),
                        Math.max(svInterval.getStart() - 500, 0),
                        Math.min(svInterval.getEnd()   + 500, sequenceDictionary.getSequence(svInterval.getContig()).getSequenceLength())))
                .map(svInterval -> makeLocatable(svInterval, sequenceDictionary))
                .collect(Collectors.toList());

        // get reads that overlap padded regions
        final JavaRDD<GATKRead> overlappingReads;
        if ( hasUserSuppliedIntervals() ) {
            overlappingReads = getReads().cache();
        } else {
            overlappingReads = myGetReads(ctx, paddedIntervals).cache();
        }

        final Broadcast<SAMSequenceDictionary> sequenceDictionaryBroadcast = ctx.broadcast(sequenceDictionary);
        final Broadcast<SVIntervalTree<Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval>>> treeBroadcast = ctx.broadcast(intervalTree);

        try {
            // pick information from each read
            final JavaPairRDD<SVInterval, AnnotatedIntervalTrio> annotatedIntervalTrios =
                    overlappingReads
                            .flatMapToPair(gatkRead -> SmallDeletionDepthCollector.collectSummaries(gatkRead, sequenceDictionaryBroadcast, treeBroadcast))
                            .reduceByKey(AnnotatedIntervalTrio::merge);

            // for now, write results as files
            SmallDeletionDepthCollector.writeBedFile(bedFile, sequenceDictionary, annotatedIntervalTrios.sortByKey().values().collect());
        } finally {
            overlappingReads.unpersist(false);
        }
    }

    // TODO: 6/8/18 this is a hack, because the interval passed in via -L will actually strip away the reads in flanking regions we want
    private JavaRDD<GATKRead> myGetReads(final JavaSparkContext ctx, final List<SimpleInterval> paddedIntervals) {

        // different way to set traversal
        final TraversalParameters traversalParameters = new TraversalParameters(paddedIntervals, false);

        final String refPath = hasReference() ? referenceArguments.getReferenceFileName() : null;

        final ReadFilter filter = makeReadFilter();

        final String readInput = readArguments.getReadFilesNames().get(0);
        return new ReadsSparkSource(ctx, readArguments.getReadValidationStringency())
                .getParallelReads(readInput, refPath, traversalParameters, bamPartitionSplitSize)
                .filter(filter::test);
    }

    // BED intervals are 0-based half open [0, ...), and GATKRead is 1-based closed [1, ...], whereas SVInterval is 0-based half-open [0, ...)
    private static SimpleInterval makeLocatable(final SVInterval svInterval, final SAMSequenceDictionary samSequenceDictionary) {
        return new SimpleInterval(samSequenceDictionary.getSequence(svInterval.getContig()).getSequenceName(), svInterval.getStart() + 1, svInterval.getEnd());
    }
}
