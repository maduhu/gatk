package org.broadinstitute.hellbender.tools.walkers.sv.filter;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

/**
 * For intervals that has a basic name annotation.
 */
@DefaultSerializer(NamedCStyleInterval.Serializer.class)
public class NamedCStyleInterval extends AnnotatedCStyleInterval {

    private final String name;

    NamedCStyleInterval(final SVInterval interval, final String name) {
        super(interval);
        this.name = name;
    }

    NamedCStyleInterval(final NamedCStyleInterval copy) {
        super(copy.interval);
        name = copy.name;
    }

    NamedCStyleInterval(final int contig, final int start, final int end, final String name) {
        this(new SVInterval(contig, start, end), name);
    }

    NamedCStyleInterval(final Kryo kryo, final Input input) {
        super(kryo, input);
        name = input.readString();
    }

    protected void serialize( final Kryo kryo, final Output output ) {
        super.serialize(kryo, output);
        output.writeString(name);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NamedCStyleInterval> {
        @Override
        public void write(final Kryo kryo, final Output output, final NamedCStyleInterval interval ) {
            interval.serialize(kryo, output);
        }

        @Override
        public NamedCStyleInterval read(final Kryo kryo, final Input input, final Class<NamedCStyleInterval> klass ) {
            return new NamedCStyleInterval(kryo, input);
        }
    }

    String getName() {
        return name;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("NamedCStyleInterval{");
        sb.append("interval=").append(interval);
        sb.append(", name='").append(name).append('\'');
        sb.append('}');
        return sb.toString();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        final NamedCStyleInterval that = (NamedCStyleInterval) o;

        return name.equals(that.name);
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + name.hashCode();
        return result;
    }
}
