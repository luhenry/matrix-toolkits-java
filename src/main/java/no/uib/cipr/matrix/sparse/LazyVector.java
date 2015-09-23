/*
 * Copyright (C) James Millard
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


	//	/**
	//	 * Upper bound on the number of potential duplicates created by <code>add()</code>.<br>
	//	 * When <code>add()</code> is called, if more unordered entries than this are present they are sorted and coalesced.
	//	 */
	//	public static final int MAX_DUPLICATES = 1 << 6;

	/**
	 * How many entries in the ordered and unordered data sections, respectively.
	 */
	int numOrdered, numUnordered;

	/**
	 * Data.
	 * 
	 * Contains <code>numOrdered</code> sorted entries at the start, and <code>numUnordered</code> entries immediately after.
	 * Unordered entries are strictly newer than ordered entries.
	 */
	double[] data;

	/**
	 * Indices to data.
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

		// TODO: should we check against zero when setting zeros?

		//TODO: ordered unordered
		int i = getIndex(index);
		data[i] = value;
	}

	@Override
	public void add(int index, double value) {
		check(index);
		//TODO: ordered unordered
		int i = getIndex(index);
		data[i] += value;
	}


	//	private int indexSearch(int index, boolean searchUnordered){
	//		int pos = Arrays.binarySearch(this.index, 0, numOrdered, index);
	//		
	//		if(pos >= 0){
	//			return pos;
	//		}
	//		
	//		pos = -1; //linear search through
	//		for(int i = numOrdered; i < numUnordered + numOrdered; i++){
	//			if(this.index[i] == index){
	//				pos = i;
	//				break;
	//			}
	//		}
	//		return pos;
	//	}

	
	/**
	 * Sorts all unordered entries, and trims array down to size. Does not remove structural zeros.
	 * Generally called as part of a method which requires sorted entries.
	 */
	private void trim(){
		if(numOrdered == index.length){return;}
		
		sort(-1);
		index = Arrays.copyOfRange(index, 0, numOrdered);
		data = Arrays.copyOfRange(data, 0, numOrdered);
	}
	
	private void sort(int newSize){
		if(numUnordered == 0) return;

		int[] tempIndex = Arrays.copyOfRange(index, numOrdered, numOrdered + numUnordered);
		double[] tempData = Arrays.copyOfRange(data, numOrdered, numOrdered + numUnordered);


		// Get a set of pointers which are ordered by the index they point to.
		Integer[] ptr = new Integer[numUnordered]; // Ugly, but, can't custom sort primitives. Java 10!
		for(int i = 0; i < ptr.length; i++){ptr[i] = Integer.valueOf(i);}
		Arrays.sort(ptr, 0, ptr.length, (ptr1, ptr2) -> {
			// negative indices in unordered area indicate an add, ignore for now.
			int ptr1Index = tempIndex[ptr1] < 0 ? ~tempIndex[ptr1] : tempIndex[ptr1];
			int ptr2Index = tempIndex[ptr2] < 0 ? ~tempIndex[ptr2] : tempIndex[ptr2];

			if(ptr1Index < ptr2Index){return -1;}
			if(ptr1Index > ptr2Index){return 1;}
			return 0; // stable sort will preserve chronological order of set()s and add()s
		});


		//Sort both data and index by index using ptrs
		for(int i = 0; i < ptr.length; i++){
			tempIndex[i] = index[ptr[i] + numOrdered];
			tempData[i] = data[ptr[i] + numOrdered];
		}

		// De-duplicate indices, process in chronological order, set()s overwrite, add()s add.
		int cursor = 0;
		for(int i = 0; i < tempIndex.length; i++){
			tempIndex[cursor] = tempIndex[cursor] < 0 ? ~tempIndex[cursor] : tempIndex[cursor];
			//TODO: de-dupe
			int j = cursor + 1;
			while( j < tempIndex.length){
				
			}
		}

		// Move to end of array to allow for reliable merge.
		if(newSize > index.length){
			int[] newIndex = new int[newSize];
			double[] newData = new double[newSize];
			System.arraycopy(index, 0, newIndex, newIndex.length - numOrdered, numOrdered);
			System.arraycopy(data, 0, newData, newData.length - numOrdered, numOrdered);
			index = newIndex;
			data = newData;
		} else {
			System.arraycopy(index, 0, index, index.length - numOrdered, numOrdered);
			System.arraycopy(data, 0, data, data.length - numOrdered, numOrdered);
		}


		merge(index, data, index.length - numOrdered, index.length,
				tempIndex, tempData, 0, cursor,
				index, data, 0);


		numOrdered = numOrdered + cursor;
		numUnordered = 0;
	}

	/**
	 * Processed from start to end. User responsible for array size checks and to ensure values
	 * cannot be overwritten if input arrays are reused for output.
	 */
	@SuppressWarnings("unused")
	private void merge(
			int[] index1, double[] data1, int start1, int end1,
			int[] index2, double[] data2, int start2, int end2,
			int[] indexOut, double[] dataOut, int startOut){


		//TODO: use starts and ends
		int i = 0, j = 0, k = 0;

		while (i < index1.length && j < index2.length){
			if (index1[i] <= index2[j]){
				indexOut[k] = index1[i];
				dataOut[k++] = data1[i++];
			}else{
				indexOut[k] = index2[j];
				dataOut[k++] = data2[j++];
			}
		}

		System.arraycopy(index1, i, indexOut, k, index1.length - i);
		System.arraycopy(data1, i, dataOut, k, data1.length - i);

		System.arraycopy(index1, j, indexOut, k, index1.length - j);
		System.arraycopy(data1, j, dataOut, k, data1.length - j);

	}

	@Override
	public double get(int index) {
		check(index);



		int in = Arrays.binarySearch(this.index, 0, used, index);
		if (in >= 0)
			return data[in];
		return 0;
	}

	/**
	 * Tries to find the index. If it is not found, a reallocation is done, and
	 * a new index is returned.
	 */
	private int getIndex(int ind) {

		int i = Arrays.binarySearch(index, 0, used, ind);

		if(i >= 0){
			return i; // Index found
		} else {
			i = ~i; // convert to insertion location
		}

		// Check available memory
		if (used >= data.length) {

			// If zero-length, use new length of 1, else double the bandwidth
			int newLength = data.length != 0 ? data.length << 1 : 1;

			// Enforce the maximum size.
			newLength = Math.min(newLength, this.size);

			// Copy existing data into new arrays
			int[] newIndex = new int[newLength];
			double[] newData = new double[newLength];

			// Copy and, make room for insertion
			System.arraycopy(index, 0, newIndex, 0, i);
			System.arraycopy(data, 0, newData, 0, i);
			System.arraycopy(index, i, newIndex, i + 1, used - i);
			System.arraycopy(data, i, newData, i + 1, used - i);

			// Update pointers
			index = newIndex;
			data = newData;
		} else if (i < used){
			// All ok, make room for insertion
			System.arraycopy(index, i, index, i + 1, used - i);
			System.arraycopy(data, i, data, i + 1, used - i);
		}


		// Insert new index
		used++;
		index[i] = ind;
		data[i] = 0.;


		// Return insertion index
		return i;
	}

	@Override
	public LazyVector copy() {
		return new LazyVector(this);
	}

	@Override
	public LazyVector zero() {
		java.util.Arrays.fill(data, 0);
		used = 0; // TODO: Confirm correctness, java doc says "preserves underlying structure" this will result it overwriting.
		return this;
	}

	@Override
	public LazyVector scale(double alpha) {
		// Quick return if possible
		if (alpha == 0)
			return zero();
		else if (alpha == 1)
			return this;

		for (int i = 0; i < used; ++i)
			data[i] *= alpha;

		return this;
	}

	@Override
	public double dot(Vector y) {
		if (!(y instanceof DenseVector))
			return super.dot(y);

		checkSize(y);

		double[] yd = ((DenseVector) y).getData();

		double ret = 0;
		for (int i = 0; i < used; ++i)
			ret += data[i] * yd[index[i]];
		return ret;
	}

	@Override
	protected double norm1() {
		double sum = 0;
		for (int i = 0; i < used; ++i)
			sum += Math.abs(data[i]);
		return sum;
	}

	@Override
	protected double norm2() {
		double norm = 0;
		for (int i = 0; i < used; ++i)
			norm += data[i] * data[i];
		return Math.sqrt(norm);
	}

	@Override
	protected double norm2_robust() {
		double scale = 0, ssq = 1;
		for (int i = 0; i < used; ++i) {
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
		double max = 0;
		for (int i = 0; i < used; ++i)
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
		return data;
	}

	/**
	 * Returns the used indices
	 */
	public int[] getIndex() {
		if (used == index.length)
			return index;

		// could run compact, or return subarray
		// compact();
		int[] indices = new int[used];
		System.arraycopy(index, 0, indices, 0, used);
		return indices;
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
		//TODO: sort
		return numOrdered;
	}

	@Override
	public void compact() {
		int nz = Matrices.cardinality(this); // catches zero entries

		if (nz < data.length) {
			int[] newIndex = new int[nz];
			double[] newData = new double[nz];

			// Copy only non-zero entries
			for (int i = 0, j = 0; i < data.length; ++i)
				if (data[i] != 0.) {
					newIndex[j] = index[i];
					newData[j] = data[i];
					j++;
				}

			data = newData;
			index = newIndex;
			used = data.length;
		}
	}

	@Override
	public Iterator<VectorEntry> iterator() {
		return new LazyVectorIterator();
	}

	@Override
	public Vector set(Vector y) {
		if (!(y instanceof SparseVector))
			return super.set(y);

		checkSize(y);

		SparseVector yc = (SparseVector) y;

		if (yc.index.length != index.length) {
			data = new double[yc.data.length];
			index = new int[yc.data.length];
		}

		System.arraycopy(yc.data, 0, data, 0, data.length);
		System.arraycopy(yc.index, 0, index, 0, index.length);
		used = yc.used;

		return this;
	}

	/**
	 * Iterator over a sparse vector
	 */
	private class LazyVectorIterator implements Iterator<VectorEntry> {

		private int cursor;

		private final LazyVectorEntry entry = new LazyVectorEntry();

		public boolean hasNext() {
			return cursor < used;
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
	 * Entry of a sparse vector
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

