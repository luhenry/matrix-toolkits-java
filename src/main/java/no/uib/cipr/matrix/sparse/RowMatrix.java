package no.uib.cipr.matrix.sparse;

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

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Iterator;
import java.util.function.Function;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.sparse.SuperIterator.SuperIteratorEntry;

/**
 * Matrix stored row-wise into vectors.
 * Arbitrary vectors can be used as rows by supplying a vector generating function which consumes <code>numColumns</code> as an int.
 */
public class RowMatrix<T extends Vector> extends AbstractMatrix {

    /**
     * Matrix data
     */
    T[] rowD;
	
	Class<T> rowClass;

    /**
     * Constructor for FlexCompRowMatrix
     * 
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of column
     */
	@SuppressWarnings("unchecked")
	public RowMatrix(int numRows, int numColumns, Function<Integer, T> supplier) {
        super(numRows, numColumns);
        
        rowClass = (Class<T>) supplier.apply(0).getClass();
        rowD =  (T[]) Array.newInstance(rowClass, numRows);
        for (int i = 0; i < numRows; ++i){
        	rowD[i] = supplier.apply(numColumns);
        }
    }
    

    /**
     * Deep Copy Constructor
     * 
     * @param A
     *            Matrix to copy contents from
     */
    @SuppressWarnings("unchecked")
	public RowMatrix(Matrix A, Function<Integer, T> supplier) {
        super(A);
        
        rowClass = (Class<T>) supplier.apply(0).getClass();
        rowD =  (T[]) Array.newInstance(rowClass, numRows);
        for (int i = 0; i < numRows; ++i){
        	rowD[i] = supplier.apply(numColumns);
        }
            
        set(A);
        
    }


    /**
     * Returns the given row
     */
    public T getRow(int i) {
        return rowD[i];
    }

    /**
     * Sets the given row equal the passed vector
     */
    public void setRow(int i, T x) {
        if (x.size() != numColumns)
            throw new IllegalArgumentException(
                    "New row must be of the same size as existing row");
        rowD[i] = x;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        checkMultAdd(x, y);

        for (int i = 0; i < numRows; ++i)
            y.add(i, alpha * rowD[i].dot(x));

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {

        checkTransMultAdd(x, y);
        
        y.scale(1. / alpha);
        
        for(int i = 0; i < rowD.length; i++){
        	y.add(x.get(i), rowD[i]);
        }
        
        return y.scale(alpha);

    }

    @Override
    public void add(int row, int column, double value) {
        rowD[row].add(column, value);
    }

    @Override
    public void set(int row, int column, double value) {
        rowD[row].set(column, value);
    }

    @Override
    public double get(int row, int column) {
        return rowD[row].get(column);
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new RowMatrixIterator();
    }

    @Override
    public Matrix copy() {
        return new FlexCompRowMatrix(this);
    }

    @Override
    public RowMatrix<T> zero() {
        for (int i = 0; i < numRows; ++i)
            rowD[i].zero();
        return this;
    }

//    @Override
//    public Matrix set(Matrix B) {} // always use element iterator

    /**
     * Tries to store the matrix as compactly as possible
     */
    public void compact() {
    	if (ISparseVector.class.isAssignableFrom(rowClass)){
            for (Vector v : rowD)
                ((ISparseVector) v).compact();
    	}
    }

    /**
     * Iterator over a matrix stored vectorwise by rows
     */
    private class RowMatrixIterator implements Iterator<MatrixEntry> {

        /**
         * Iterates over each row vector
         */
        private SuperIterator<T, VectorEntry> iterator = new SuperIterator<T, VectorEntry>(
                Arrays.asList(rowD));

        
        
        /**
         * Entry returned
         */
        private RowMatrixEntry entry = new RowMatrixEntry();

        public boolean hasNext() {
        	
            return iterator.hasNext();
        }

        public MatrixEntry next() {
            SuperIteratorEntry<VectorEntry> se = iterator.next();
            entry.update(se.index(), se.get());
            return entry;
        }

        public void remove() {
            iterator.remove();
        }

    }

    /**
     * Entry of a matrix stored vectorwise by rows
     */
    private static class RowMatrixEntry implements MatrixEntry {

        private int row;

        private VectorEntry entry;

        public void update(int row, VectorEntry entry) {
            this.row = row;
            this.entry = entry;
        }

        public int row() {
            return row;
        }

        public int column() {
            return entry.index();
        }

        public double get() {
            return entry.get();
        }

        public void set(double value) {
            entry.set(value);
        }

    }

}

