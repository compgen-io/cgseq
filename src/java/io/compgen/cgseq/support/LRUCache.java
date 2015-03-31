package io.compgen.cgseq.support;

import java.util.LinkedHashMap;
import java.util.Map.Entry;

public class LRUCache<K,V> {
	private class LRUMap extends LinkedHashMap <K,V> {
		private static final long serialVersionUID = -6904223367386557660L;
		final private int maxSize;
		
		public LRUMap(int maxSize) {
			this.maxSize = maxSize;
		}

		@Override
		protected boolean removeEldestEntry(Entry<K, V> eldest) {
			return (size() > maxSize);
		}
	}

	final private LRUMap map;
	public LRUCache(int maxSize) {
		map = new LRUMap(maxSize);
	}
	
	public boolean containsKey(K key) {
		return map.containsKey(key);
	}
	
	public V get(K key) {
		if (map.containsKey(key)) {
			V value = map.remove(key);
			map.put(key, value);
			return value;
		}
		return null;
	}
	
	public void put(K key, V value) {
		map.put(key, value);
	}
}