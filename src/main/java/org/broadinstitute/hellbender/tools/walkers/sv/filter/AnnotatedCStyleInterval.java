package org.broadinstitute.hellbender.tools.walkers.sv.filter;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

/**
 * Hope SVInterval is not final, but I understand it's trying to be lightweight....
 *
 * For use by adding various annotations as necessary.
 *
 * Interval must be C-style 0-based semi-open interval (e.g. SVInterval, BED interval),
 * as apposed to {@link htsjdk.samtools.util.Locatable}'s, which are [1, ...].
 */
public abstract class AnnotatedCStyleInterval {
    private static final SVInterval.Serializer intervalSerializer = new SVInterval.Serializer();

    protected final SVInterval interval;

    public AnnotatedCStyleInterval(final SVInterval interval) {
        this.interval = interval;
    }

    protected AnnotatedCStyleInterval(final Kryo kryo, final Input input) {
        interval = intervalSerializer.read(kryo, input, SVInterval.class);
    }

    protected void serialize( final Kryo kryo, final Output output ) {
        intervalSerializer.write(kryo, output, interval);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final AnnotatedCStyleInterval that = (AnnotatedCStyleInterval) o;

        return interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        return interval.hashCode();
    }

    // delegation ======================================================================================================

    public int getContig() {
        return interval.getContig();
    }

    public int getStart() {
        return interval.getStart();
    }

    public int getEnd() {
        return interval.getEnd();
    }

    public int getLength() {
        return interval.getLength();
    }

    /**
     * Check {@link SVInterval#overlaps(SVInterval)}
     */
    public boolean overlaps( final SVInterval that ) {
        return interval.overlaps(that);
    }

    /**
     * See {@link #overlaps(SVInterval)}.
     */
    public boolean overlaps( final AnnotatedCStyleInterval that ) {
        return interval.overlaps(that.interval);
    }
}
