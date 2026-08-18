// spectrum.go extends the kmer package with spectral analysis utilities:
// k-mer spectrum computation, frequency filtering, and Jaccard similarity
// between k-mer sets.
package kmer

import (
	"math"
	"sort"
)

// Spectrum represents the frequency distribution of k-mers.
type Spectrum struct {
	K       int
	Counts  map[string]int
	Total   int // total number of k-mers observed
	Unique  int // number of distinct k-mers
}

// NewSpectrum computes the k-mer spectrum for a sequence.
func NewSpectrum(seq string, k int) (*Spectrum, error) {
	counts, err := Count(seq, k)
	if err != nil {
		return nil, err
	}
	total := 0
	for _, c := range counts {
		total += c
	}
	return &Spectrum{K: k, Counts: counts, Total: total, Unique: len(counts)}, nil
}

// TopN returns the n most frequent k-mers sorted by count descending.
func (s *Spectrum) TopN(n int) []KmerCount {
	sorted := s.Sorted()
	if n > len(sorted) {
		n = len(sorted)
	}
	return sorted[:n]
}

// Sorted returns all k-mers sorted by count descending, then by string.
func (s *Spectrum) Sorted() []KmerCount {
	out := make([]KmerCount, 0, len(s.Counts))
	for km, c := range s.Counts {
		out = append(out, KmerCount{Kmer: km, Count: c})
	}
	sort.Slice(out, func(a, b int) bool {
		if out[a].Count != out[b].Count {
			return out[a].Count > out[b].Count
		}
		return out[a].Kmer < out[b].Kmer
	})
	return out
}

// KmerCount pairs a k-mer string with its frequency.
type KmerCount struct {
	Kmer  string
	Count int
}

// FilterAbove returns k-mers with count >= threshold.
func (s *Spectrum) FilterAbove(threshold int) map[string]int {
	out := make(map[string]int)
	for km, c := range s.Counts {
		if c >= threshold {
			out[km] = c
		}
	}
	return out
}

// FilterBelow returns k-mers with count < threshold.
func (s *Spectrum) FilterBelow(threshold int) map[string]int {
	out := make(map[string]int)
	for km, c := range s.Counts {
		if c < threshold {
			out[km] = c
		}
	}
	return out
}

// Frequency returns the relative frequency of a specific k-mer.
func (s *Spectrum) Frequency(kmer string) float64 {
	if s.Total == 0 {
		return 0
	}
	return float64(s.Counts[kmer]) / float64(s.Total)
}

// Entropy computes the Shannon entropy of the k-mer distribution.
func (s *Spectrum) Entropy() float64 {
	if s.Total == 0 {
		return 0
	}
	var h float64
	for _, c := range s.Counts {
		p := float64(c) / float64(s.Total)
		if p > 0 {
			h -= p * math.Log2(p)
		}
	}
	return h
}

// JaccardSimilarity computes the Jaccard index between two k-mer sets.
// J(A,B) = |A ∩ B| / |A ∪ B|
func JaccardSimilarity(a, b map[string]int) float64 {
	if len(a) == 0 && len(b) == 0 {
		return 1.0
	}
	intersect := 0
	for km := range a {
		if _, ok := b[km]; ok {
			intersect++
		}
	}
	union := len(a)
	for km := range b {
		if _, ok := a[km]; !ok {
			union++
		}
	}
	if union == 0 {
		return 0
	}
	return float64(intersect) / float64(union)
}

// CosineSimilarity computes the cosine similarity between two k-mer count
// vectors.
func CosineSimilarity(a, b map[string]int) float64 {
	dot := 0.0
	magA := 0.0
	magB := 0.0
	all := make(map[string]bool)
	for km := range a {
		all[km] = true
	}
	for km := range b {
		all[km] = true
	}
	for km := range all {
		ca := float64(a[km])
		cb := float64(b[km])
		dot += ca * cb
		magA += ca * ca
		magB += cb * cb
	}
	if magA == 0 || magB == 0 {
		return 0
	}
	return dot / (math.Sqrt(magA) * math.Sqrt(magB))
}
