package kmer

import (
	"math"
	"testing"
)

func TestNewSpectrum(t *testing.T) {
	s, err := NewSpectrum("ACGTACGT", 2)
	if err != nil {
		t.Fatalf("spectrum: %v", err)
	}
	if s.K != 2 {
		t.Fatalf("K = %d", s.K)
	}
	if s.Total != 7 {
		t.Fatalf("Total = %d, want 7", s.Total)
	}
	if s.Unique == 0 {
		t.Fatal("Unique should be > 0")
	}
}

func TestSpectrumTopN(t *testing.T) {
	s, _ := NewSpectrum("AAAA", 2)
	top := s.TopN(1)
	if len(top) != 1 || top[0].Kmer != "AA" || top[0].Count != 3 {
		t.Fatalf("top = %v", top)
	}
}

func TestSpectrumFilter(t *testing.T) {
	s, _ := NewSpectrum("ACGTACGT", 2)
	above := s.FilterAbove(2)
	below := s.FilterBelow(2)
	totalAbove := 0
	for range above {
		totalAbove++
	}
	totalBelow := 0
	for range below {
		totalBelow++
	}
	if totalAbove+totalBelow != s.Unique {
		t.Fatalf("above(%d) + below(%d) != unique(%d)", totalAbove, totalBelow, s.Unique)
	}
}

func TestSpectrumEntropy(t *testing.T) {
	s, _ := NewSpectrum("AAAA", 1)
	// All A's: entropy should be 0 (only one k-mer).
	if s.Entropy() != 0 {
		t.Fatalf("entropy = %f, want 0", s.Entropy())
	}
	s2, _ := NewSpectrum("ACGT", 1)
	// Equal distribution: entropy should be 2 bits.
	if math.Abs(s2.Entropy()-2.0) > 0.01 {
		t.Fatalf("entropy = %f, want 2.0", s2.Entropy())
	}
}

func TestJaccardSimilarity(t *testing.T) {
	a := map[string]int{"AC": 1, "CG": 1, "GT": 1}
	b := map[string]int{"AC": 2, "CG": 1, "TT": 1}
	j := JaccardSimilarity(a, b)
	// intersection {AC,CG} = 2, union {AC,CG,GT,TT} = 4 → 0.5
	if math.Abs(j-0.5) > 0.01 {
		t.Fatalf("jaccard = %f, want 0.5", j)
	}
}

func TestCosineSimilarity(t *testing.T) {
	a := map[string]int{"AC": 2, "CG": 0}
	b := map[string]int{"AC": 2, "CG": 0}
	c := CosineSimilarity(a, b)
	if math.Abs(c-1.0) > 0.01 {
		t.Fatalf("cosine = %f, want 1.0", c)
	}
}
