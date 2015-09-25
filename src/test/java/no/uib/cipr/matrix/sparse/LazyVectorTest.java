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

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Utilities;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.VectorTestAbstract;

/**
 * Test of LazyVector
 */
public class LazyVectorTest extends VectorTestAbstract {

    public LazyVectorTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Math.min(Utilities.getInt(max), n);
        x = new LazyVector(n);
        xd = Utilities.populate(x, m);
    }

    public void testLazyVectorIndices() {
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

    public void testBug27() {
        double[] tfVector = {0.0, 0.5, 0.0, 0.4, 0.0};
        DenseVector dense = new DenseVector(tfVector, false);
        LazyVector vectorTF = new LazyVector(dense);
        vectorTF.compact();

        int used = vectorTF.getUsed();
        
        assertTrue(used == 2); // vectorTF.getUsed() returns 5

        for (Iterator<VectorEntry> it = vectorTF.iterator(); it.hasNext();) {
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
    public void testOverAllocation() {
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

    public void testGetRawIndex() {
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

    public void testGetRawData() {
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
    
    public void testPerformance(){
    	testPerformance(1000);
    	testPerformance(1000);
    	testPerformance(1000);
    	System.out.println();
    	System.out.println("Fill\tSparse\tLazy\tHash\tTree");
    	
    	for(double fill = 2<<14; fill >= 1; fill = fill/2){
    		testPerformance((int)fill);
    	}
    }
    
    
    public void testPerformance(int fill){
    	
    	int size = fill*10;
    	int mult = 1;
    	int loops = 100000;
    	
    	while(loops*(long)(fill) > 1E7){
    		loops /= 2;
    		mult *= 2;
    	}
    	

    	
    	double[][] data = new double[loops][fill];
    	int[][] index = new int[loops][fill];
    	for(int i = 0; i < loops; i++){
    		for(int f = 0; f < fill; f++){
    			data[i][f] = Math.random();
    			index[i][f] = (int)(Math.random()*size);
    		}
    	}
    	
//    	int warmUp = Math.min(100000, loops);
//    	for(int i = 0; i < warmUp; i++){
//    		LazyVector lv = new LazyVector(size);
//    		SparseVector sv = new SparseVector(size);
//    		MapVector hv = new MapVector(size, HashMap::new);
//    		MapVector tv = new MapVector(size, TreeMap::new);
//    		for(int f = 0; f < fill; f++){
//    			lv.set(index[i%loops][f], data[i%loops][f]);
//    			sv.set(index[i%loops][f], data[i%loops][f]);
//    			hv.set(index[i%loops][f], data[i%loops][f]);
//    			tv.set(index[i%loops][f], data[i%loops][f]);
//    		}
//    	}


    	long lazyStart = System.nanoTime();
    	for(int i = 0; i < loops; i++){
    		LazyVector lv = new LazyVector(size);
    		for(int f = 0; f < fill; f++){
    			lv.set(index[i][f], data[i][f]);
    		}
    		//lv.compact();
    	}
    	long lazyEnd = System.nanoTime();
    	
    	long sparseStart = System.nanoTime();
    	for(int i = 0; i < loops; i++){
    		SparseVector sv = new SparseVector(size);
    		for(int f = 0; f < fill; f++){
    			sv.set(index[i][f], data[i][f]);
    		}
    		//sv.compact();
    	}
    	long sparseEnd = System.nanoTime();
    	
    	long hashStart = System.nanoTime();
    	for(int i = 0; i < loops; i++){
    		MapVector hv = new MapVector(size, HashMap::new);
    		for(int f = 0; f < fill; f++){
    			hv.set(index[i][f], data[i][f]);
    		}
    		//hv.compact();
    	}
    	long hashEnd = System.nanoTime();

    	
    	long treeStart = System.nanoTime();
    	for(int i = 0; i < loops; i++){
    		MapVector tv = new MapVector(size, TreeMap::new);
    		for(int f = 0; f < fill; f++){
    			tv.set(index[i][f], data[i][f]);
    		}
    		//tv.compact();
    	}
    	long treeEnd = System.nanoTime();
    	
//    	System.out.printf("%d\t", fill);
//    	System.out.printf("Sparse\t%f\t", mult*(sparseEnd - sparseStart)/1000000.0);
//    	System.out.printf("Lazy\t%f\t", mult*(lazyEnd - lazyStart)/1000000.0);
//    	System.out.printf("Hash\t%f\t", mult*(hashEnd - hashStart)/1000000.0);
//    	System.out.printf("Tree\t%f\t", mult*(treeEnd - treeStart)/1000000.0);
    	System.out.printf("%d\t%f\t%f\t%f\t%f\t", 
    			fill, 
    			mult*(sparseEnd - sparseStart)/1000000.0, 
    			mult*(lazyEnd - lazyStart)/1000000.0,
    			mult*(hashEnd - hashStart)/1000000.0,
    			mult*(treeEnd - treeStart)/1000000.0);
    	System.out.println();
    }
    
}
