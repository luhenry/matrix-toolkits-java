/**
 * May 21, 2007
 * @author Samuel Halliday, ThinkTank Maths Limited
 * Copyright ThinkTank Maths Limited 2007
 */
package no.uib.cipr.matrix.sparse;

import no.uib.cipr.matrix.Vector;

/**
 * @author Samuel Halliday, ThinkTank Maths Limited
 */
public interface ISparseVector extends Vector {

    /**
     * Returns the indices
     */
    public int[] getIndex();

    /**
     * Number of entries used in the sparse structure
     */
    public int getUsed();
    
    /**
     * Returns the internal value array. This array may contain extra elements
     * beyond the number that are used. If it is greater than the number used,
     * the remaining values will be 0. Since this vector can resize its internal
     * data, if it is modified, this array may no longer represent the internal
     * state.
     * 
     * @return The internal array of values.
     */
    public double[] getData();
    
    /**
     * Compacts the vector
     */
    public void compact();
}
