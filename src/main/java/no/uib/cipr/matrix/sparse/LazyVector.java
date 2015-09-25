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
	 * If true, allow adds to be put in the unordered area with bitwise inverse indices.
	 * If false an <code>add()</code> that cannot be combined within <code>COALAESE_DEPTH</code>
	 * will trigger a sort.
	 */
	public static final boolean ALLOW_DUPLICATE_UNORDERED_ADDS = true;
	
	/**
	 * Depth to search the most recent end of the unordered section
	 * in an attempt to combine adds and overwrite with sets.
	 */
	// Approximately 1-2 cache lines
	public static final int COALESCE_DEPTH = 16;

	
	/**
	 * Max number of entries in unordered section that <code>get()</code> will perform a linear search on.
	 * If get is called and unordered is larger than this all entries will be converted to ordered.
	 */
	// must be <= COALESCE_DEPTH, otherwise it can't be assumed that unordered entries are unique
	// must be >= 0.
	public static final int MAX_GET_DEPTH = 16; 
	
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


		int in = Arrays.binarySearch(this.index, 0, numOrdered, index);
		if (in >= 0){
			return data[in];
		}
			
		if(numUnordered > MAX_GET_DEPTH){
			convertAllToOrdered();
			return get(index);
		}
		
		
		in = partialLinearSearch(index, MAX_GET_DEPTH);
		
		if (in >= 0){
			return data[in];
		}
		
		return 0;
	}

	
	/**
	 * @param ind
	 * @param val
	 * @param doAdd
	 */
	private void put(int ind, double val, boolean doAdd){ 
		// Search ordered, if found here it is guaranteed not to be in unordered section.
		//int in = Arrays.binarySearch(index, 0, numOrdered, ind);

		int in = no.uib.cipr.matrix.sparse.Arrays.binarySearchSmaller(index, ind);
		
		
		//Search recent end of unordered.
		in = (in >= 0 && index[in] != ind)? in : partialLinearSearch(ind, COALESCE_DEPTH);
		
		
		if(in >= 0){
			if(doAdd){
				data[in] += val; //combine with sets or adds leaving index or inverse index in place.
			}else{
				data[in] = val;
				index[in] = ind; // if in tail of unordered, adds are overwritten as sets.
			}
			return;
		}
		
		
		// if ordered > 3/4, reallocate here. After this block there is guaranteed to be non-zero room above the ordered section
		if(index.length == 0 || numOrdered > (index.length/2)){
			int newSize = Math.min(size, Math.max(8, index.length*2));
			int[] newIndex = new int[newSize];
			double[] newData = new double[newSize];
			System.arraycopy(index, 0, newIndex, 0, numOrdered + numUnordered);
			System.arraycopy(data, 0, newData, 0, numOrdered + numUnordered);
			index = newIndex;
			data = newData;
		}
		
		// If no more room in unordered area compact and retry/reallocate.
		// remove second condition to allow duplicate adds in unordered area
		if(numOrdered + numUnordered == index.length || doAdd && numUnordered > COALESCE_DEPTH){ 
			convertAllToOrdered();
			put(ind, val, doAdd);
			return;
		}
		
		// Append to ordered section if possible
		if(numUnordered == 0 && (numOrdered == 0 || ind > index[numOrdered - 1])){
			data[numOrdered] = val;
			index[numOrdered] = ind;
			numOrdered++;
			return;
		}

		
		data[numOrdered + numUnordered] = val;
		index[numOrdered + numUnordered] = doAdd ? ~ind : ind;		
		numUnordered++;
	}
	
	
	/**
	 * Reverse linear search through the unordered section. Matches on index, or the bitwise inverse of index.
	 * 
	 * @return The position of the last matching entry in the index array, or -1 if not found.
	 */
	private int partialLinearSearch(int ind, int depth){
		int inv = ~ind;
		
		int start = numUnordered + numOrdered - 1;
		int end = Math.max(start - depth, numOrdered - 1);
		
		for(int i = start; i > end; i--){
			if(index[i] == ind || index[i] == inv){
				return i;
			}
		}
		
		return -1;
	}
	
	/**
	 * Sorts all unordered entries, and trims array down to size. Does not remove structural zeros.
	 * Generally called as part of a method which requires sorted entries.
	 */
	private void trim(){
		if(numOrdered == index.length){return;}
		
		convertAllToOrdered();
		index = Arrays.copyOfRange(index, 0, numOrdered);
		data = Arrays.copyOfRange(data, 0, numOrdered);
	}
	
	/** Sorts and coalesces unordered entries, and combines
	 * @param newSize used to dictate size of new array if resize is necessary. If less than 0 a new array will not be allocated.
	 */
	private void convertAllToOrdered(){
		if(numUnordered == 0) return;
		
		//sort but preserve chronological order of sets and adds
		int[] tempIndex = new int[numUnordered];
		double[] tempData = new double[numUnordered];
		copyAndSortUnordered(tempIndex, tempData);

		// coalesce sets and adds
		int cursor = fullCoalesce(tempIndex, tempData);


		// Move to end of array to allow for reliable merge, allocate new arrays if requested.
		System.arraycopy(index, 0, index, index.length - numOrdered, numOrdered);
		System.arraycopy(data, 0, data, data.length - numOrdered, numOrdered);


		// Merge old ordered entries and newly sorted/coalesced "unordered" entries.
		numOrdered = merge(index, data, index.length - numOrdered, index.length,
				tempIndex, tempData, 0, cursor,
				index, data, 0);

		numUnordered = 0;

		Arrays.fill(index, numOrdered, index.length, 0); // not really needed. Could just let junk data sit there.
	}

	/** 
	 * Accepts unordered entries and performs a stable sort by index.
	 * Bitwise inverse indices are sorted as though they are positive.
	 * Requires that argument arrays are copies of the 
	 */
	private void copyAndSortUnordered(int[] tempIndex, double[] tempData) {
		// Get a set of pointers which are ordered by the index they point to.
		
		Integer[] ptr = new Integer[numUnordered]; // Ugly, but, can't custom sort primitives. Java 10!
		
		for(int i = 0; i < numUnordered; i++){ptr[i] = Integer.valueOf(i + numOrdered);}
		
		Arrays.sort(ptr, 0, ptr.length, (ptr1, ptr2) -> {
			// negative indices in unordered area indicate an add, ignore for now.
			int ptr1Index = index[ptr1] < 0 ? ~index[ptr1] : index[ptr1];
			int ptr2Index = index[ptr2] < 0 ? ~index[ptr2] : index[ptr2];

			if(ptr1Index < ptr2Index){return -1;}
			if(ptr1Index > ptr2Index){return 1;}
			return 0; // stable sort will preserve chronological order of set()s and add()s
		});


		//Sort both data and index by index using ptrs
		for(int i = 0; i < ptr.length; i++){
			tempIndex[i] = index[ptr[i]];
			tempData[i] = data[ptr[i]];
		}
	}

	
	/** Performs a full compaction of the history stored in the unordered entries
	 * Input must be sorted by index, and for equal indices sorted chronologically
	 * 
	 * @return new length of data after coalesce of adds and sets
	 */
	private int fullCoalesce(int[] index, double[] data){
		int cursor1 = 0;
		int cursor2 = 0;
		
		while(cursor2 < index.length){
			
			index[cursor1] = index[cursor2] < 0 ? ~index[cursor2] : index[cursor2]; // start new run, if first is an "add" convert to set.
			data[cursor1] =  data[cursor2];
			cursor2++;
			
			int inv = ~index[cursor1];
			
			while(cursor2 < index.length && (index[cursor2] == index[cursor1] || index[cursor2] == inv)){
				
				if(index[cursor2] < 0){
					data[cursor1] += data[cursor2];
				} else {
					data[cursor1] = data[cursor2];
				}
					
				cursor2++;
			}
			
			cursor1++;
		}
				
		
		return cursor1;
	}
	
	
	/**Takes two pairs of index sorted array slices and does an index sorted merge.
	 * Processed from start to end, entries duplicated in both inputs are removed with preference to keep entries from index2 and data2.
	 * User responsible for array size checks and to ensure values
	 * cannot be overwritten if input arrays are reused for output.
	 * 
	 * @return length of data after
	 */
	private int merge(
			int[] index1, double[] data1, int start1, int end1,
			int[] index2, double[] data2, int start2, int end2,
			int[] indexOut, double[] dataOut, int startOut){

	
		int i = start1, j = start2, k = startOut;

		while (i < end1 && j < end2){
			if (index1[i] < index2[j]){
				indexOut[k] = index1[i];
				dataOut[k++] = data1[i++];
			}else if (index1[i] > index2[j]){
				indexOut[k] = index2[j];
				dataOut[k++] = data2[j++];
			} else {
				indexOut[k] = index2[j];
				dataOut[k++] = data2[j++];
				i++;
				//throw new IllegalArgumentException("Index present in both index inputs. Inputs must be unique.");
			}
		}

		System.arraycopy(index1, i, indexOut, k, end1 - i);
		System.arraycopy(data1, i, dataOut, k, end1 - i);

		System.arraycopy(index2, j, indexOut, k, end2 - j);
		System.arraycopy(data2, j, dataOut, k, end2 - j);

		
		return (k- startOut) + (end1 - i) + (end2 - j);
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
		
		convertAllToOrdered();
		// could run compact, or return subarray
		// compact();
		int[] indices = new int[numOrdered];
		System.arraycopy(index, 0, indices, 0, numOrdered);
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
		convertAllToOrdered();
		return numOrdered;
	}

	@Override
	public void compact() {
		
		
		//TODO: loop through with two cursors and ditch zeros, then call trim;
		
		convertAllToOrdered();
		int cursor = 0;
		for(int i = 0; i < numOrdered; i++){
			if(data[i] != 0){
				data[cursor] = data[i];
				index[cursor] = index[i];
				cursor++;
			}
		}
		numOrdered = cursor;
		
		trim();
		
	}

	@Override
	public Iterator<VectorEntry> iterator() {
		convertAllToOrdered();
		return new LazyVectorIterator();
	}

	@Override
	public Vector set(Vector y) {
//		if (!(y instanceof SparseVector))
			return super.set(y);

//		checkSize(y);
//
//		SparseVector yc = (SparseVector) y;
//
//		if (yc.index.length != index.length) {
//			data = new double[yc.data.length];
//			index = new int[yc.data.length];
//		}
//
//		System.arraycopy(yc.data, 0, data, 0, data.length);
//		System.arraycopy(yc.index, 0, index, 0, index.length);
//		used = yc.used;
//
//		return this;
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

