/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import com.google.common.collect.Ordering;

import java.util.Objects;

/**
 * Based on response in http://stackoverflow.com/questions/2670982/using-pairs-or-2-tuples-in-java
 * @author Eli Rodgers-Melnick
 */
public class Tuple<X, Y> implements Comparable<Tuple<X, Y>> {
    public final X x;
    public final Y y;
    /**
     * Instantiates a tuple object, which just holds 2 values
     * @param x The first object
     * @param y The second object
     */
    public Tuple(X x, Y y) {
        this.x = x;
        this.y = y;
    }

    public X getX() {
        return x;
    }

    public Y getY() {
        return y;
    }

    @Override
    public int hashCode() {
        return (x.hashCode()^y.hashCode());
    }
    @Override
    public String toString() {
        return "("+x.toString()+","+y.toString()+")";
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Tuple<X, Y> other = (Tuple<X, Y>) obj;
        if (!Objects.equals(this.x, other.x)) {
            return false;
        }
        if (!Objects.equals(this.y, other.y)) {
            return false;
        }
        return true;
    }

    @Override
    public int compareTo(Tuple<X, Y> o) {
        if(x instanceof Double) return Double.compare((Double)x, (Double)o.x);

        if(x instanceof Comparable) {
            return ((Comparable) x).compareTo(o.x);
        }
        return 0;
    }
}
