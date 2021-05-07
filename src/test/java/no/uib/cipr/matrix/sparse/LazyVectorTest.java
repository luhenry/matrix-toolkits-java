/*
 * Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 * 
 * This file is part of MTJ.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package no.uib.cipr.matrix.sparse;

import no.uib.cipr.matrix.*;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.function.BiFunction;
import java.util.function.Function;

import static org.junit.Assert.*;

/**
 * Test of LazyVector
 */
@Ignore
public class LazyVectorTest extends VectorTestAbstract {

    public LazyVectorTest() {
        super();
    }

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Math.min(Utilities.getInt(max), n);
        x = new LazyVector(n);
        xd = Utilities.populate(x, m);
    }

    @Test
    @Ignore
    public void testLazyVectorIndices() throws Exception {
        /*
         * MTJ subtlety in getIndex() for SparseVector. before calling
         * getIndex(), you must call compact()... implementations may choose to
         * do nothing in this call, but the Intel extended LAPACK
         * implementations (and MTJ's SparseVector) require it. An alternative
         * to vector.getIndex() is VectorMethods.getIndex(Vector) which will
         * wrap this for you. It can take an arbitrary Vector and if it can be
         * cast to a SparseVector will compact it and use its getIndex() method
         * instead. (just so you're aware of this). Sam.
         */

        // check that "infinite dimensions" doesn't use infinite memory
        LazyVector vector = new LazyVector(Integer.MAX_VALUE);
        int[] index = vector.getIndex();
        assert index != null;
        assert index.length == 0;

        // check that creating with double[] with zeros works
        double[] entries = new double[5];
        entries[0] = 0.0;
        entries[1] = 0.0;
        entries[2] = 1.0;
        entries[3] = 0.0;
        entries[4] = 2.0;
        Vector dense = new DenseVector(entries, false);
        vector = new LazyVector(dense);

        // NOTE: must compact before calling getIndex()!!!
        // vector.compact();
        index = vector.getIndex();
        assert index != null;
        assert index.length == 5 : "expected length of 5, but got "
                + index.length + ", with elements " + Arrays.toString(index);
    }

    @Test
    public void testBug27()  throws Exception {
        double[] tfVector = {0.0, 0.5, 0.0, 0.4, 0.0};
        DenseVector dense = new DenseVector(tfVector, false);
        LazyVector vectorTF = new LazyVector(dense);
        vectorTF.compact();

        int used = vectorTF.getUsed();

        assertTrue(used == 2); // vectorTF.getUsed() returns 5

        for (Iterator<VectorEntry> it = vectorTF.iterator(); it.hasNext(); ) {
            VectorEntry ve = it.next();
            int index = ve.index();
            double value = ve.get();
            assertTrue(tfVector[index] == value);
        }
    }

    /**
     * Unit test checking that the LazyVector does not end up ever using more
     * than "size" elements.
     */
    @Test
    public void testOverAllocation() throws Exception {
        for (int d = 0; d < 10; d++) {
            LazyVector v = new LazyVector(d, 0);
            assertEquals(0, v.index.length);
            assertEquals(0, v.data.length);

            // Fill with non-zero elements.
            for (int i = 0; i < d; i++) {
                v.set(i, 1.0 + i);
            }

            assertEquals(d, v.index.length);
            assertEquals(d, v.data.length);
        }
    }

    @Test
    public void testGetRawIndex() throws Exception {
        LazyVector vector = new LazyVector(Integer.MAX_VALUE);
        int[] index = vector.getRawIndex();
        assertTrue(index != null);
        assertTrue(index.length == 0);
        assertSame(index, vector.index);
        assertEquals(index.length, vector.getRawData().length);

        vector.set(2, 1.0);
        vector.set(1, 0.0);
        vector.set(4, 2.0);

        index = vector.getRawIndex();
        assertSame(index, vector.index);
        assertEquals(index.length, vector.getRawData().length);

        // In this case, the raw index is larger than the used, so the raw
        // indices have more entries than the other one.
        assertTrue(index.length >= vector.getUsed());
        assertTrue(index.length >= vector.getIndex().length);
    }

    @Test
    public void testGetRawData() throws Exception {
        LazyVector vector = new LazyVector(Integer.MAX_VALUE);
        double[] data = vector.getRawData();
        assertTrue(data != null);
        assertTrue(data.length == 0);
        assertSame(data, vector.data);
        assertEquals(data.length, vector.getRawIndex().length);

        vector.set(2, 1.0);
        vector.set(1, 0.0);
        vector.set(4, 2.0);

        data = vector.getRawData();
        assertSame(data, vector.data);
        assertEquals(data.length, vector.getRawIndex().length);

        // In this case, the raw index is larger than the used, so the raw
        // indices have more entries than the other one.
        assertTrue(data.length >= vector.getUsed());
        assertTrue(data.length >= vector.getIndex().length);
    }

    @Test
    @Ignore
    public void testPerformance()  throws Exception {
        testPerformance(100, 10, 5);


        System.out.println();
        System.out.println("Fill\tSparse\tLazy\tHash\tTree");

        for (double fill = 2 << 8; fill >= 1; fill = fill / 2) {
            testPerformance((int) fill, 10, 5);
        }
    }


    private void testPerformance(int fill, int sizeFact, int repeats) {

        int size = fill * sizeFact;
        double mult = 0.1;
        int loops = 1000000;

        while (loops * (long) (fill) > 1E7) {
            loops /= 2;
            mult *= 2;
        }

        long sparse = Long.MAX_VALUE, lazy = Long.MAX_VALUE, hash = Long.MAX_VALUE, tree = Long.MAX_VALUE;

        for (int repeat = 0; repeat < repeats; repeat++) {
            double[][] data = new double[loops][fill];
            int[][] index = new int[loops][fill];
            for (int i = 0; i < loops; i++) {
                for (int f = 0; f < fill; f++) {
                    data[i][f] = Math.random();
                    index[i][f] = (int) (Math.random() * size);
                }
            }

            sparse = Math.min(sparse, trial(index, data, size, SparseVector::new));

            lazy = Math.min(lazy, trial(index, data, size, LazyVector::new));
            hash = Math.min(hash, trial(index, data, size, dim -> new MapVector(dim, HashMap::new)));
            tree = Math.min(tree, trial(index, data, size, dim -> new MapVector(dim, TreeMap::new)));
        }

        System.out.printf("%d\t%f\t%f\t%f\t%f\t",
                fill,
                mult * (sparse) / 1000000.0,
                mult * (lazy) / 1000000.0,
                mult * (hash) / 1000000.0,
                mult * (tree) / 1000000.0);
        System.out.println();
    }

    private long trial(int[][] index, double[][] data, int size, Function<Integer, ISparseVector> supplier) {
        System.gc();
        long start = System.nanoTime();
        for (int i = 0; i < index.length; i++) {
            ISparseVector v = supplier.apply(size);
            for (int f = 0; f < index[i].length; f++) {
                v.set(index[i][f], data[i][f]);
            }
            //v.compact();
            //v.getData();
        }
        return System.nanoTime() - start;
    }

    private long trialPreAllocate(int[][] index, double[][] data, int size, BiFunction<Integer, Integer, ISparseVector> supplier) {
        long start = System.nanoTime();
        for (int i = 0; i < index.length; i++) {
            ISparseVector v = supplier.apply(size, index[i].length);
            for (int f = 0; f < index[i].length; f++) {
                v.set(index[i][f], data[i][f]);
            }
            //v.compact();
            //v.getData();
        }

        return System.nanoTime() - start;
    }
}
