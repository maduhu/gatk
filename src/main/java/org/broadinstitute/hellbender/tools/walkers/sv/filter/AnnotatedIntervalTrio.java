package org.broadinstitute.hellbender.tools.walkers.sv.filter;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

/**
 * A trio of contiguous small intervals annotated with "depth" information,
 * currently designed specifically for filtering small (< ~ 1kbp) deletions.
 *
 * By "contiguous" we mean the left interval's end base neighbors middle interval's start base,
 * and similarly for middle and right intervals.
 */
@DefaultSerializer(AnnotatedIntervalTrio.Serializer.class)
public final class AnnotatedIntervalTrio {
    private static final DepthAnnotatedCStyleInterval.Serializer serializer = new DepthAnnotatedCStyleInterval.Serializer();
    private final DepthAnnotatedCStyleInterval leftFlank;
    private final DepthAnnotatedCStyleInterval middle;
    private final DepthAnnotatedCStyleInterval rightFlank;

    AnnotatedIntervalTrio(final DepthAnnotatedCStyleInterval leftFlank, final DepthAnnotatedCStyleInterval middle,
                          final DepthAnnotatedCStyleInterval rightFlank) {
        this.leftFlank = leftFlank;
        this.middle = middle;
        this.rightFlank = rightFlank;
    }

    AnnotatedIntervalTrio(final GATKRead gatkRead, final SVInterval readInterval,
                          final Tuple3<NamedCStyleInterval, NamedCStyleInterval, NamedCStyleInterval> leftMiddleAndRight) {
        this(   new DepthAnnotatedCStyleInterval(leftMiddleAndRight._1()),
                new DepthAnnotatedCStyleInterval(leftMiddleAndRight._2()),
                new DepthAnnotatedCStyleInterval(leftMiddleAndRight._3()));

        leftFlank.collectFromRead(gatkRead, readInterval);
        middle.collectFromRead(gatkRead, readInterval);
        rightFlank.collectFromRead(gatkRead, readInterval);
    }

    AnnotatedIntervalTrio(final Kryo kryo, final Input input) {
        this.leftFlank = kryo.readObject(input, DepthAnnotatedCStyleInterval.class, serializer);
        this.middle = kryo.readObject(input, DepthAnnotatedCStyleInterval.class, serializer);
        this.rightFlank = kryo.readObject(input, DepthAnnotatedCStyleInterval.class, serializer);
    }

    SVInterval getFlankedCStyleInterval() {
        return new SVInterval(middle.getContig(), leftFlank.getStart(), rightFlank.getEnd());
    }

    String toBed_5_String(final SAMSequenceDictionary sequenceDictionary) {
        return leftFlank.toBed_5_String(sequenceDictionary) + middle.toBed_5_String(sequenceDictionary) + rightFlank.toBed_5_String(sequenceDictionary);
    }

    String toSummaryString(final SAMSequenceDictionary sequenceDictionary) {
        return leftFlank.toSummaryString(sequenceDictionary) + middle.toSummaryString(sequenceDictionary) + rightFlank.toSummaryString(sequenceDictionary);
    }

    public static AnnotatedIntervalTrio merge(final AnnotatedIntervalTrio first, final AnnotatedIntervalTrio second) {
        return new AnnotatedIntervalTrio(
                DepthAnnotatedCStyleInterval.merge(first.leftFlank, second.leftFlank),
                DepthAnnotatedCStyleInterval.merge(first.middle, second.middle),
                DepthAnnotatedCStyleInterval.merge(first.rightFlank, second.rightFlank));
    }

    private void serialize( final Kryo kryo, final Output output ) {
        kryo.writeObject(output, leftFlank,  serializer);
        kryo.writeObject(output, middle,     serializer);
        kryo.writeObject(output, rightFlank, serializer);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AnnotatedIntervalTrio> {
        @Override
        public void write( final Kryo kryo, final Output output, final AnnotatedIntervalTrio trio ) {
            trio.serialize(kryo, output);
        }

        @Override
        public AnnotatedIntervalTrio read(final Kryo kryo, final Input input, final Class<AnnotatedIntervalTrio> klass ) {
            return new AnnotatedIntervalTrio(kryo, input);
        }
    }
}
