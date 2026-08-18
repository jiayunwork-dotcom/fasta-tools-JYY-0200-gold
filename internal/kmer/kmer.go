package kmer

import (
	"fmt"
)

// Count returns overlapping k-mer frequencies of s.
// k must be >= 1 and len(s) must be >= k. An empty or too-short sequence
// yields a non-nil empty map.
func Count(s string, k int) (map[string]int, error) {
	if k < 1 {
		return nil, fmt.Errorf("k must be >= 1, got %d", k)
	}
	if len(s) < k {
		return map[string]int{}, nil
	}
	counts := make(map[string]int)
	for i := 0; i+k <= len(s); i++ {
		counts[s[i:i+k]]++
	}
	return counts, nil
}
