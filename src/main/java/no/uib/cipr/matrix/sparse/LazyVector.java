/*
 * Based on SparseVector Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 * 
 * Copyright (C) 2015 James Millard 
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

import java.util.Iterator;
import java.util.Arrays;

import no.uib.cipr.matrix.AbstractVector;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;

/**
 * A sparse vector. Allows sets and adds to be performed in constant time by delaying sorting until required.
 * Internal arrays sort and resize when data is requested.
 * Alternating often between accessing and inserting data will result in poor performance.
 */
@SuppressWarnings("serial")
public class LazyVector extends AbstractVector implements ISparseVector {


	/**
	 * How many entries in the ordered and unordered data sections, respectively.
	 */
	int numOrdered, numUnordered;

	/**
	 * Data. <br>
	 * Contains <code>numOrdered</code> sorted entries at the start, and <code>numUnordered</code> entries immediately after.
	 * Unordered entries are strictly newer than ordered entries.
	 */
	double[] data;

	/**
	 * Indices to data. <br>
	 * Contains <code>numOrdered</code> sorted entries at the start, and <code>numUnordered</code> entries immediately after.
	 * Unordered entries are strictly newer than ordered entries.
	 * In the unordered section, <code>add()</code>s are stored using the bitwise inverse of their index,
	 * resulting in strictly negative values.
	 */
	int[] index;



	/**
	 * Constructor for SparseVector.
	 * 
	 * @param size
	 *            Size of the vector
	 * @param nz
	 *            Initial number of non-zeros
	 */
	public LazyVector(int size, int nz) {
		super(size);
		data = new double[nz];
		index = new int[nz];
	}

	/**
	 * Constructor for SparseVector. Zero initial pre-allocation
	 * 
	 * @param size
	 *            Size of the vector
	 */
	public LazyVector(int size) {
		this(size, 0);
	}

	/**
	 * Constructor for SparseVector, and copies the contents from the supplied
	 * vector. Zero initial pre-allocation
	 * 
	 * @param x
	 *            Vector to copy from. A deep copy is made
	 */
	public LazyVector(Vector x) {
		super(x);

		int nz = Matrices.cardinality(x);
		data = new double[nz];
		index = new int[nz];
		set(x);
	}



	/**
	 * Constructor for SparseVector
	 * 
	 * @param size
	 *            Size of the vector
	 * @param index
	 *            Indices of the vector
	 * @param data
	 *            Entries of the vector
	 * @param deep
	 *            True for a deep copy. For shallow copies, the given indices
	 *            will be used internally
	 */
	public LazyVector(int size, int[] index, double[] data, boolean deep) {
		super(size);

		if (index.length != data.length)
			throw new IllegalArgumentException("index.length != data.length");

		if (deep) {
			this.index = index.clone();
			this.data = data.clone();
		} else {
			this.index = index;
			this.data = data;
		}

		boolean ordered = true;
		for(int i = 1; i < index.length; i++){
			if(index[i] < index[i-1]){
				ordered = false;
				break;
			}
		}

		if(ordered){
			numOrdered = index.length;
		} else {
			numUnordered = index.length;
		}
	}

	/**
	 * Constructor for SparseVector
	 * 
	 * @param size
	 *            Size of the vector
	 * @param index
	 *            The vector indices are copies from this array
	 * @param data
	 *            The vector entries are copies from this array
	 */
	public LazyVector(int size, int[] index, double[] data) {
		this(size, index, data, true);
	}

	@Override
	public void set(int index, double value) {
		check(index);
		put(index, value, false);
	}

	@Override
	public void add(int index, double value) {
		check(index);
		put(index, value, true);
	}


	@Override
	public double get(int index) {
		check(index);

		convertAllToOrdered();
		int in = Arrays.binarySearch(this.index, 0, numOrdered, index);
		return in < 0 ? 0 : data[in];
	}
	
		
	/**
	 * @param ind
	 * @param val
	 * @param doAdd
	 */
	private void put(int ind, double val, boolean doAdd){ 

		if(numUnordered == 0){

			int pos = ~numOrdered;
			if(numOrdered > 0 && ind <= index[numOrdered - 1]){ // If cannot add to end
				pos = Arrays.binarySearch(index, 0, numOrdered, ind);
				
				if(pos >= 0){
					data[pos] = doAdd ? data[pos] + val : val;
					return;
				}
			}


			if(val == 0){ return;}

			// This entire insertion block is optional, it just improves performance for small N.
			pos = ~pos; //bitwise inverse indicates insertion is necessary, convert to insertion position
			if(numOrdered - pos < 128){

				int[] newIndex = index;
				double[] newData = data;

				// Check available memory
				if (numOrdered >= data.length) {

					// If zero-length, use new length of 1, else double the bandwidth
					int newLength = Math.min(size, Math.max(4, index.length*2));

					// Copy existing data into new arrays
					newIndex = new int[newLength];
					newData = new double[newLength];
					System.arraycopy(index, 0, newIndex, 0, pos);
					System.arraycopy(data, 0, newData, 0, pos);
				}

				// All ok, make room for insertion
				System.arraycopy(index, pos, newIndex, pos + 1, numOrdered - pos);
				System.arraycopy(data, pos, newData, pos + 1, numOrdered - pos);

				// Put in new structure
				newIndex[pos] = ind;
				newData[pos] = val;
				numOrdered++;

				index = newIndex;
				data = newData;
				return;
			}

		}


		// If no more room in unordered area compact and retry. Reallocate if numOrdered is high.
		if(numOrdered + numUnordered == index.length ){

			int newSize = -1;
			if(numOrdered > (index.length/2) || index.length == 0 ){
				newSize = Math.min(size, Math.max(4,   index.length*2));
			}
			convertAllToOrdered(newSize);
			put(ind, val, doAdd);
			return;
		}


		data[numOrdered + numUnordered] = val;
		index[numOrdered + numUnordered] = doAdd? ~ind : ind;
		numUnordered++;
	}



	private void convertAllToOrdered(){
		convertAllToOrdered(-1);
	}

	/**
	 * Combines all unordered entries into the ordered section.
	 * @param newSize used to dictate size of new array if resize is necessary. If less than 0 a new array will not be allocated.
	 */
	private void convertAllToOrdered(int newSize){
		if(numUnordered == 0 && index.length >= newSize) return;

		newSize = Math.max(newSize, index.length);

		// Put a sorted copy of unordered entries at the end of the new array
		int[] newIndex = new int[newSize];
		double[] newData = new double[newSize];
		copyAndSortUnordered(
				numOrdered, index, data,
				newIndex.length - numUnordered, newIndex, newData,
				numUnordered);

		// Merge old ordered entries and coalesce and merge new unordered entries.
		numOrdered = merge(index, data, 0, numOrdered,
				newIndex, newData, newIndex.length - numUnordered, newIndex.length,
				newIndex, newData, 0);

		numUnordered = 0;

		index = newIndex;
		data = newData;
	}

	/** 
	 * Transfers a sorted compy of the unordered entries into the argument arrays.
	 * Bitwise inverse indices (adds) are sorted as though they are their positive equivalent.
	 * Requires that argument arrays are large enough to hold all unordered entries
	 */

	private static void copyAndSortUnordered(
			int startIn, int[] indexIn, double[] dataIn,
			int startOut, int[] indexOut, double[] dataOut,
			int length) {

		Integer[] ptr = new Integer[length]; // Ugly, but, can't custom sort primitives. Java 10?!
		int[] positiveIndex = new int[length]; //Arrays.copyOfRange(index, numOrdered, numOrdered + numUnordered);
		for(int i = 0; i < length; i++){
			ptr[i] = i;
			positiveIndex[i] = Math.max(indexIn[i+startIn], ~indexIn[i+startIn]);
		}

		Arrays.sort(ptr, 0, ptr.length, (ptr1, ptr2) -> {
			int d = positiveIndex[ptr1] - positiveIndex[ptr2]; // Can't underflow, both positive
			return d;
			//return d == 0 ? ptr1.compareTo(ptr2) : d;
		});


		//Sort both data and index by index using ptrs
		for(int i = 0; i < ptr.length; i++){
			indexOut[startOut + i] = indexIn[ptr[i]+startIn];
			dataOut[startOut + i] = dataIn[ptr[i]+startIn];
		}

		return;
	}

	/**Takes two pairs of index sorted array slices and does an index sorted merge.
	 * Processed from start to end, entries duplicated in both inputs are removed with preference to keep entries from index2 and data2.
	 * User responsible for array size checks and to ensure values
	 * cannot be overwritten if input arrays are reused for output. TODO update doc to include coalesce functionality
	 * 
	 * @return length of data after
	 */
	private static int merge(
			int[] index1, double[] data1, int start1, int end1,
			int[] index2, double[] data2, int start2, int end2,
			int[] indexOut, double[] dataOut, int startOut){


		int i = start1, j = start2, k = startOut;

		while (i < end1 && j < end2){
			int ind1 = index1[i];
			int ind2 = index2[j] < 0 ? ~index2[j] : index2[j];

			if (ind1 < ind2){
				indexOut[k] = ind1;
				dataOut[k] = data1[i++];
			}else if (ind1 > ind2){
				indexOut[k] = ind2;
				dataOut[k] = data2[j++];
				j = collectDuplicateUnorderedEntries(index2, data2, end2, dataOut, j, k, ind2); 
			} else {
				indexOut[k] = ind1;
				dataOut[k] = data1[i++];
				j = collectDuplicateUnorderedEntries(index2, data2, end2, dataOut, j, k, ind2);
			}

			if(dataOut[k] != 0){k++;} // removes 0s from ordered entries
		}

		while(i < end1){
			indexOut[k] = index1[i];
			dataOut[k] = data1[i++];
			if(dataOut[k] != 0){k++;}
		}

		while(j < end2){
			int ind2 = index2[j] < 0 ? ~index2[j] : index2[j];
			indexOut[k] = ind2;
			dataOut[k] = data2[j++];
			j = collectDuplicateUnorderedEntries(index2, data2, end2, dataOut, j, k, ind2);
			if(dataOut[k] != 0){k++;}
		}

		return k - startOut;
	}

	/**
	 * Sub-method of merge.
	 * If there are unordered entries with the same index collect and apply them to dataOut.
	 * @return "j", final position in unordered entries array
	 */
	private static int collectDuplicateUnorderedEntries(int[] index2, double[] data2, int end2, double[] dataOut, int j, int k, int ind2) {
		int inv = ~ind2;
		while(j < end2 && (index2[j] == ind2 || index2[j] == inv)){
			dataOut[k] = index2[j] < 0 ? dataOut[k] + data2[j] : data2[j];
			j++;
		}
		return j;
	}

	@Override
	public LazyVector copy() {
		return new LazyVector(this);
	}

	@Override
	public LazyVector zero() {
		java.util.Arrays.fill(data, 0);
		numOrdered = 0; // TODO: Confirm correctness, java doc says "preserves underlying structure" this will result it overwriting.
		numUnordered = 0;
		return this;
	}

	@Override
	public LazyVector scale(double alpha) {
		// Quick return if possible
		if (alpha == 0)
			return zero();
		else if (alpha == 1)
			return this;

		convertAllToOrdered();
		for (int i = 0; i < numOrdered; ++i)
			data[i] *= alpha;

		return this;
	}

	@Override
	public double dot(Vector y) {
		if (!(y instanceof DenseVector))
			return super.dot(y);

		checkSize(y);
		convertAllToOrdered();

		double[] yd = ((DenseVector) y).getData();

		double ret = 0;
		for (int i = 0; i < numOrdered; ++i) //TODO: can unroll to encourage better use of SSE/AVX
			ret += data[i] * yd[index[i]];
		return ret;
	}

	@Override
	protected double norm1() {
		convertAllToOrdered();
		double sum = 0;
		for (int i = 0; i < numOrdered; ++i)
			sum += Math.abs(data[i]);
		return sum;
	}

	@Override
	protected double norm2() {
		convertAllToOrdered();
		double norm = 0;
		for (int i = 0; i < numOrdered; ++i)
			norm += data[i] * data[i];
		return Math.sqrt(norm);
	}

	@Override
	protected double norm2_robust() {
		convertAllToOrdered();
		double scale = 0, ssq = 1;
		for (int i = 0; i < numOrdered; ++i) {
			if (data[i] != 0) {
				double absxi = Math.abs(data[i]);
				if (scale < absxi) {
					ssq = 1 + ssq * Math.pow(scale / absxi, 2);
					scale = absxi;
				} else
					ssq = ssq + Math.pow(absxi / scale, 2);
			}
		}
		return scale * Math.sqrt(ssq);
	}

	@Override
	protected double normInf() {
		convertAllToOrdered();
		double max = 0;
		for (int i = 0; i < numOrdered; ++i)
			max = Math.max(Math.abs(data[i]), max);
		return max;
	}

	/**
	 * Returns the internal value array. This array may contain extra elements
	 * beyond the number that are used. If it is greater than the number used,
	 * the remaining values will be 0. Since this vector can resize its internal
	 * data, if it is modified, this array may no longer represent the internal
	 * state.
	 * 
	 * @return The internal array of values.
	 */
	public double[] getData() {
		convertAllToOrdered();
		return data;
	}

	/**
	 * Returns the used indices
	 */
	public int[] getIndex() {
		if (numOrdered == index.length)
			return index;

		compact();
		return index;
		//		
		//		convertAllToOrdered();
		//		// could run compact, or return subarray
		//		
		//		int[] indices = new int[numOrdered];
		//		System.arraycopy(index, 0, indices, 0, numOrdered);
		//		return indices;
	}

	/**
	 * Gets the raw internal index array. This array may contain extra elements
	 * beyond the number that are used. If it is greater than the number used,
	 * the remaining indices will be 0. Since this vector can resize its
	 * internal data, if it is modified, this array may no longer represent the
	 * internal state.
	 * 
	 * @return The internal array of indices, whose length is greater than or
	 *         equal to the number of used elements. Indices in the array beyond
	 *         the used elements are not valid indices since they are unused.
	 */
	public int[] getRawIndex() {
		return index;
	}

	/**
	 * Gets the raw internal data array. This array may contain extra elements
	 * beyond the number that are used. If it is greater than the number used,
	 * the remaining indices will be 0. Since this vector can resize its
	 * internal data, if it is modified, this array may no longer represent the
	 * internal state.
	 * 
	 * @return The internal array of values, whose length is greater than or
	 *         equal to the number of used elements. Values in the array beyond
	 *         the used elements are not valid since they are unused.
	 */
	public double[] getRawData() {
		return data;
	}

	/**
	 * Number of entries used in the sparse structure
	 */
	public int getUsed() {
		convertAllToOrdered();
		return numOrdered;
	}

	@Override
	public void compact() {

		convertAllToOrdered();
		if(numOrdered == index.length){return;}

		int cursor = 0;
		for(int i = 0; i < numOrdered; i++){
			if(data[i] != 0){
				data[cursor] = data[i];
				index[cursor] = index[i];
				cursor++;
			}
		}
		numOrdered = cursor;

		index = Arrays.copyOfRange(index, 0, numOrdered);
		data = Arrays.copyOfRange(data, 0, numOrdered);
	}

	@Override
	public Iterator<VectorEntry> iterator() {
		convertAllToOrdered();
		return new LazyVectorIterator();
	}

	@Override
	public Vector set(Vector y) {
		if (!(y instanceof ISparseVector)){
			return super.set(y);
		}


		checkSize(y);

		ISparseVector yc = (ISparseVector) y;

		int[] ycInd = yc.getIndex();
		double[] ycData = yc.getData();


		if (ycInd.length != index.length) {
			data = new double[ycData.length];
			index = new int[ycData.length];
		}

		System.arraycopy(ycData, 0, data, 0, data.length);
		System.arraycopy(ycInd, 0, index, 0, index.length);


		numUnordered = yc.getUsed(); // If getIndex is strictly ordered then this can be changed

		return this;
	}

	/**
	 * Iterator over a LazyVector
	 */
	private class LazyVectorIterator implements Iterator<VectorEntry> {

		private int cursor;

		private final LazyVectorEntry entry = new LazyVectorEntry();

		public boolean hasNext() {
			return cursor < numOrdered;
		}

		public VectorEntry next() {
			entry.update(cursor);

			cursor++;

			return entry;
		}

		public void remove() {
			entry.set(0);
		}

	}

	/**
	 * Entry of a LazyVector
	 */
	private class LazyVectorEntry implements VectorEntry {

		private int cursor;

		public void update(int cursor) {
			this.cursor = cursor;
		}

		public int index() {
			return index[cursor];
		}

		public double get() {
			return data[cursor];
		}

		public void set(double value) {
			data[cursor] = value;
		}

	}

}

