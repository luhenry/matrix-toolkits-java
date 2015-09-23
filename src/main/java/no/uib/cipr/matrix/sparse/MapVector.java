package no.uib.cipr.matrix.sparse;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.function.Supplier;
import java.util.Arrays;

import no.uib.cipr.matrix.AbstractVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;


/**
 * Sparse Vector implemented using standard library HashMap
 */
public class MapVector extends AbstractVector implements ISparseVector{

	final Map<Integer, Double> map;
	Supplier<Map<Integer, Double>> supplier;
	
    public MapVector(int size) {
        this(size, HashMap::new);
    }
	
    public MapVector(int size, Supplier<Map<Integer, Double>> supplier) {
        super(size);
        this.supplier = supplier;
        map = supplier.get();
    }
    
    public MapVector(Vector x){
    	this(x, HashMap::new);
    }
    
    public MapVector(Vector x, Supplier<Map<Integer, Double>> supplier) {
        super(x.size());
        this.supplier = supplier;
        map = supplier.get();
        x.forEach(e ->{
        	map.put(e.index(), e.get());
        }); 
    }
	    
	public int size() {
		return map.size();
	}

	@Override
	public void set(int index, double value) {
		check(index);
		map.put(index, value);
	}

	@Override
	public void add(int index, double value) {
		check(index);
		double current = map.getOrDefault(index, 0.0);
		map.put(index, value + current);
	}

	@Override
	public double get(int index) {
		check(index);
		return map.getOrDefault(index, 0.0);
	}

	@Override
	public Vector copy() {
		return new MapVector(this, supplier);
	}

	@Override
	public Vector zero() {
//		for(Entry<Integer, Double> e: map.entrySet()){
//			e.setValue(0.0);
//		}
		map.clear();
		return this;
	}

	@Override
	public Vector scale(double alpha) {
		map.keySet().forEach(k ->{
			map.put(k, map.get(k)*alpha);
		});
		return this;
	}

	@Override
	public Vector set(Vector y) {
		return set(1, y);
	}

	@Override
	public Vector set(double alpha, Vector y) {
		map.clear();
		y.forEach(e ->{
			map.put(e.index(), e.get()*alpha);
		});
		return null;
	}

	@Override
	public Vector add(Vector y) {
		return add(1, y);
	}

	@Override
	public Vector add(double alpha, Vector y) {
		y.forEach(e ->{
			map.put(e.index(), map.getOrDefault(e.index(), 0.0) + e.get()*alpha);
		});
		return this;
	}

	@Override
	public double dot(Vector y) {
		double sum = 0;
		map.entrySet().stream().map(e -> e.getValue()*y.get(e.getKey()));
		return sum;
	}

    @Override
    public Iterator<VectorEntry> iterator() {
        return new MapVectorIterator();
    }

    /**
     * Iterator over a sparse vector
     */
    private class MapVectorIterator implements Iterator<VectorEntry> {

        private int cursor;
        private int[] indices = MapVector.this.getIndex();
        private final MapVectorEntry entry = new MapVectorEntry();

        public boolean hasNext() {
            return cursor < indices.length;
        }

        public VectorEntry next() {
            entry.update(indices[cursor]);

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
    private class MapVectorEntry implements VectorEntry {

        private int index;

        public void update(int index) {
            this.index = index;
        }

        public int index() {
            return index;
        }

        public double get() {
            return map.getOrDefault(index, 0.0);
        }

        public void set(double value) {
            map.put(index, value);
        }

    }

	@Override
	public int[] getIndex() {
		int[] index = new int[map.size()];
		int i = 0;
		
		for(Entry<Integer, Double> e: map.entrySet()){
			index[i++] = e.getKey();
		}
		Arrays.sort(index);
		
		return index;
	}

	@Override
	public int getUsed() {
		return map.size();
	}

	@Override
	public void compact() {
		
		Map<Integer, Double> newMap = supplier.get();
		
		map.entrySet().stream()
		.filter(e -> e.getValue() != 0)
		.forEach(e -> newMap.put(e.getKey(), e.getValue()));
		
		map.clear();
		map.putAll(newMap);

	}

}
