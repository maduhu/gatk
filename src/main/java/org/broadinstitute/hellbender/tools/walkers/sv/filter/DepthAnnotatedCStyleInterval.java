package org.broadinstitute.hellbender.tools.walkers.sv.filter;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A special type of annotated interval.
 *
 * Note that if the interval is large (heuristically speaking, over 1~2 kbp),
 * one should consider if CNV bins is a better option.
 */
@DefaultSerializer(DepthAnnotatedCStyleInterval.Serializer.class)
final class DepthAnnotatedCStyleInterval extends NamedCStyleInterval {
    private final int[] readStartCount; // count of reads that starts at this index into {@code interval} (expected to be mostly 1 after mark duplicate)
    private final int[] baseCount; // count of bases at this index into {@code interval}
    private long totalMQ;          // accumulation of MQ of short reads that starts within {@code interval}

    private DepthAnnotatedCStyleInterval(final SVInterval interval, final String name) {
        this(new NamedCStyleInterval(interval, name));
    }

    DepthAnnotatedCStyleInterval(final NamedCStyleInterval interval) {
        super(interval);

        this.readStartCount = new int[interval.getLength()];
        this.baseCount = new int[interval.getLength()];
        this.totalMQ = 0;
    }

    DepthAnnotatedCStyleInterval(final Kryo kryo, final Input input) {
        super(kryo, input);

        int length = input.readInt();
        readStartCount = new int[length];
        for (int i = 0; i < length; ++i) {readStartCount[i] = input.readInt();}
        length = input.readInt();
        baseCount = new int[length];
        for (int i = 0; i < length; ++i) {baseCount[i] = input.readInt();}
        totalMQ = input.readLong();
    }

    protected void serialize( final Kryo kryo, final Output output ) {
        super.serialize(kryo, output);

        output.writeInt(readStartCount.length);
        for (final int cnt : readStartCount) {
            output.writeInt(cnt);
        }
        output.writeInt(baseCount.length);
        for (final int cnt : baseCount) {
            output.writeInt(cnt);
        }
        output.writeLong(totalMQ);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<DepthAnnotatedCStyleInterval> {
        @Override
        public void write( final Kryo kryo, final Output output, final DepthAnnotatedCStyleInterval interval ) {
            interval.serialize(kryo, output);
        }

        @Override
        public DepthAnnotatedCStyleInterval read(final Kryo kryo, final Input input, final Class<DepthAnnotatedCStyleInterval> klass ) {
            return new DepthAnnotatedCStyleInterval(kryo, input);
        }
    }

    // core ============================================================================================================

    void collectFromRead(final GATKRead gatkRead, final SVInterval readInterval) {
        bumpBaseCount(readInterval);
        bumpReadStartCount(readInterval);
        addReadMQ(readInterval, gatkRead.getMappingQuality());
    }

    private void bumpBaseCount(final SVInterval readInterval) {
        if (interval.overlaps(readInterval)) {
            final int start = Math.max(0, readInterval.getStart() - interval.getStart());
            final int end = Math.min(readInterval.getEnd(), interval.getEnd()) - interval.getStart();
            for (int i = start; i < end; ++i) {
                ++baseCount[i];
            }
        }
    }

    private void bumpReadStartCount(final SVInterval readInterval) {
        if ( interval.overlaps(readInterval) && readInterval.getStart() >= interval.getStart() )
            ++ readStartCount[ readInterval.getStart() - interval.getStart() ];
    }

    private void addReadMQ(final SVInterval readInterval, final int mq) {
        if ( interval.overlaps(readInterval) && readInterval.getStart() >= interval.getStart() )
            totalMQ += mq;
    }

    static DepthAnnotatedCStyleInterval merge(final DepthAnnotatedCStyleInterval one, final DepthAnnotatedCStyleInterval two) {
        Utils.validate(one.interval.equals(two.interval),
                "trying to merge two annotated intervals attached to different intervals.\tone: " +
                        one.interval.toString() + " two: " + two.interval.toString());
        final DepthAnnotatedCStyleInterval merged = new DepthAnnotatedCStyleInterval(one.interval, one.getName());
        for (int i = 0; i < one.interval.getLength(); ++i) {
            merged.readStartCount[i] = one.readStartCount[i] + two.readStartCount[i];
            merged.baseCount[i] = one.baseCount[i] + two.baseCount[i];
            merged.totalMQ = one.totalMQ + two.totalMQ;
        }
        return merged;
    }

    // output ==========================================================================================================

    private int getMeanMQ() {
        int totalReadCount = 0;
        for (final int cnt : readStartCount) {
            totalReadCount += cnt;
        }
        return totalReadCount == 0 ? 0 : (int) totalMQ/totalReadCount;
    }

    /**
     * 4th column, i.e. the NAME column, contains: [meanMQ;readCnt;baseCnt;name]
     * 5th column contains the mean MQ.
     */
    String toBed_5_String(final SAMSequenceDictionary samSequenceDictionary) {
        return String.format("%s\t%d\t%d\t%s\t%d\n",
                samSequenceDictionary.getSequence(interval.getContig()).getSequenceName(),
                interval.getStart(),
                interval.getEnd(),
                "Q" + getMeanMQ() + ";" + "R" + MathUtils.sum(readStartCount) + ";" + "B" + MathUtils.median(baseCount) + ";" + getName(),
                interval.getLength()) // score
                ;
    }

    // BED-like string, where the 4th and 5th column are some arrays (so won't work nice with IGV, but not intended to either)
    String toSummaryString(final SAMSequenceDictionary samSequenceDictionary) {
        return String.format("%s\t%d\t%d\t%s\t%s\n",
                samSequenceDictionary.getSequence(interval.getContig()).getSequenceName(),
                interval.getStart(),
                interval.getEnd(),
                joinIntArrayWithComma(readStartCount),
                joinIntArrayWithComma(baseCount))
                ;
    }

    private static String joinIntArrayWithComma(final int[] array) {
        StringBuilder stringBuilder = new StringBuilder(2 * array.length);
        stringBuilder.append(array[0]);
        for (int i = 1; i < array.length; ++i) {stringBuilder.append(",").append(array[i]);}
        return stringBuilder.toString();
    }

}
