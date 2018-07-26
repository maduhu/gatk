package org.broadinstitute.hellbender.tools.walkers.sv.filter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import htsjdk.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "TO BE WRITTEN",
        summary = "TO BE WRITTEN",
        programGroup = VariantFilteringProgramGroup.class)
public class CollectDepthAnnotationForSmallIntervals extends FeatureWalker<SimpleBEDFeature> {

    @Argument(doc = "path to input BED file holding intervals of interest (trailing columns will be kept)", shortName = "bed", fullName = "inputBED")
    private String inputBED;

    @Argument(doc = "path to write BED file for the intervals and flanking regions",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String bedFile;

    private SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
    private final List<AnnotatedIntervalTrio> annotatedIntervalTrioList = new ArrayList<>(7_000); // capacity is a guess

    //==================================================================================================================

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
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(SimpleBEDFeature.class);
    }

    @Override
    public Path getDrivingFeatureFile() {
        return IOUtils.getPath(inputBED);
    }

    /**
     * Traverse the input simple BED features, use reads that overlap each
     * (flanked via {@link SmallDeletionDepthCollector#getEqualSizeFlankedInterval(SAMSequenceDictionary, SimpleBEDFeature)})
     * interval from the feature, to add annotation.
     */
    @Override
    public void traverse() {
        CountingReadFilter readFilter = makeReadFilter();
        // Process each feature in the input stream.
        Utils.stream(drivingFeatures).forEach(feature -> {

            final SimpleInterval equalSizeFlankedInterval = SmallDeletionDepthCollector.getEqualSizeFlankedInterval(sequenceDictionary, feature);

            final ReadsContext overlappingReads = new ReadsContext(getReads(), equalSizeFlankedInterval, readFilter);

            apply(feature,
                    overlappingReads,
                    new ReferenceContext(getReference(), equalSizeFlankedInterval),
                    getFeatureContext(equalSizeFlankedInterval));
            progressMeter.update(feature);
        });
    }

    // use reads for each feature
    @Override
    public void apply(final SimpleBEDFeature feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        annotatedIntervalTrioList.add(SmallDeletionDepthCollector.forWalkerApply(feature, readsContext, referenceContext, featureContext, sequenceDictionary));
    }

    @Override
    public Object onTraversalSuccess() {
        annotatedIntervalTrioList.sort(Comparator.comparing(AnnotatedIntervalTrio::getFlankedCStyleInterval));
        SmallDeletionDepthCollector.writeBedFile(bedFile, sequenceDictionary, annotatedIntervalTrioList);
        return super.onTraversalSuccess();
    }
}
