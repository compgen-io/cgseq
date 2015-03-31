package io.compgen.cgseq.support;

import io.compgen.common.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MapCount<T extends Comparable<T>> {
	private Map<T, Integer> counter = new HashMap<T, Integer>();
	
	public void incr(T k) {
		if (counter.containsKey(k)) {
			counter.put(k, counter.get(k) + 1);
		} else {
			counter.put(k, 1);
		}
	}
	
	public List<Pair<T, Integer>> getSortedCounts() {
		List<Pair<T, Integer>> out = new ArrayList<Pair<T, Integer>>();
		
		for (T k: counter.keySet()) {
			out.add(new Pair<T, Integer>(k, counter.get(k)));
		}
		
		Collections.sort(out, new Comparator<Pair<T, Integer>>() {

			@Override
			public int compare(Pair<T, Integer> o1,
					Pair<T, Integer> o2) {
				// reverse order
				return -Integer.compare(o1.two, o2.two);
			}});
		
		return out;
	}
}
