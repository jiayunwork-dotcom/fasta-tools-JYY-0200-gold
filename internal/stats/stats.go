// Package stats provides statistical analysis utilities for FASTA sequence
// collections. It computes sequence length distributions, base composition,
// N50/L50 metrics, and other assembly-quality statistics commonly used in
// genomics and bioinformatics.
package stats

import (
	"errors"
	"math"
	"sort"

	"fasta-tools/internal/fasta"
)

// ErrNoRecords is returned when the input has no sequences.
var ErrNoRecords = errors.New("stats: no records")

// Summary holds computed statistics for a set of FASTA records.
type Summary struct {
	NumRecords    int
	TotalBases    int
	MinLength     int
	MaxLength     int
	MeanLength    float64
	MedianLength  float64
	N50           int
	L50           int
	GCPercent     float64
	BaseFreq      BaseFrequency
}

// BaseFrequency holds per-base counts.
type BaseFrequency struct {
	A, C, G, T, U, N, Other int
}

// Total returns the total number of bases counted.
func (bf *BaseFrequency) Total() int {
	return bf.A + bf.C + bf.G + bf.T + bf.U + bf.N + bf.Other
}

// Compute calculates statistics from a set of FASTA records.
func Compute(records []fasta.Record) (*Summary, error) {
	if len(records) == 0 {
		return nil, ErrNoRecords
	}
	s := &Summary{NumRecords: len(records)}
	lengths := make([]int, len(records))
	for i, r := range records {
		ln := len(r.Sequence)
		lengths[i] = ln
		s.TotalBases += ln
		countBases(r.Sequence, &s.BaseFreq)
	}
	sort.Ints(lengths)
	s.MinLength = lengths[0]
	s.MaxLength = lengths[len(lengths)-1]
	s.MeanLength = float64(s.TotalBases) / float64(s.NumRecords)
	s.MedianLength = median(lengths)
	s.N50, s.L50 = computeN50L50(lengths, s.TotalBases)
	gc := s.BaseFreq.G + s.BaseFreq.C
	total := s.BaseFreq.Total()
	if total > 0 {
		s.GCPercent = float64(gc) / float64(total) * 100
	}
	return s, nil
}

// countBases updates the frequency table from a sequence string.
func countBases(seq string, bf *BaseFrequency) {
	for i := 0; i < len(seq); i++ {
		switch seq[i] {
		case 'A', 'a':
			bf.A++
		case 'C', 'c':
			bf.C++
		case 'G', 'g':
			bf.G++
		case 'T', 't':
			bf.T++
		case 'U', 'u':
			bf.U++
		case 'N', 'n':
			bf.N++
		default:
			bf.Other++
		}
	}
}

// median computes the median of a sorted slice.
func median(sorted []int) float64 {
	n := len(sorted)
	if n == 0 {
		return 0
	}
	if n%2 == 0 {
		return float64(sorted[n/2-1]+sorted[n/2]) / 2.0
	}
	return float64(sorted[n/2])
}

// computeN50L50 calculates N50 and L50 from sorted lengths.
func computeN50L50(sorted []int, totalBases int) (int, int) {
	half := totalBases / 2
	cumulative := 0
	// Traverse from longest to shortest.
	for i := len(sorted) - 1; i >= 0; i-- {
		cumulative += sorted[i]
		if cumulative >= half {
			return sorted[i], len(sorted) - i
		}
	}
	return 0, 0
}

// LengthHistogram builds a histogram of sequence lengths with the given bin
// width.
func LengthHistogram(records []fasta.Record, binWidth int) []HistogramBin {
	if len(records) == 0 || binWidth <= 0 {
		return nil
	}
	lengths := make([]int, len(records))
	maxLen := 0
	for i, r := range records {
		lengths[i] = len(r.Sequence)
		if lengths[i] > maxLen {
			maxLen = lengths[i]
		}
	}
	numBins := maxLen/binWidth + 1
	bins := make([]int, numBins)
	for _, l := range lengths {
		bin := l / binWidth
		bins[bin]++
	}
	var out []HistogramBin
	for i, count := range bins {
		if count > 0 {
			out = append(out, HistogramBin{
				Lower: i * binWidth,
				Upper: (i + 1) * binWidth,
				Count: count,
			})
		}
	}
	return out
}

// HistogramBin represents a single bin in a length histogram.
type HistogramBin struct {
	Lower int
	Upper int
	Count int
}

// Percentile returns the p-th percentile of sequence lengths (0 <= p <= 100).
func Percentile(records []fasta.Record, p float64) float64 {
	if len(records) == 0 || p < 0 || p > 100 {
		return 0
	}
	lengths := make([]int, len(records))
	for i, r := range records {
		lengths[i] = len(r.Sequence)
	}
	sort.Ints(lengths)
	idx := (p / 100) * float64(len(lengths)-1)
	lo := int(math.Floor(idx))
	hi := int(math.Ceil(idx))
	if lo == hi || hi >= len(lengths) {
		return float64(lengths[lo])
	}
	frac := idx - float64(lo)
	return float64(lengths[lo])*(1-frac) + float64(lengths[hi])*frac
}

// FilterByLength returns only records whose length falls within [min, max].
func FilterByLength(records []fasta.Record, min, max int) []fasta.Record {
	var out []fasta.Record
	for _, r := range records {
		l := len(r.Sequence)
		if l >= min && l <= max {
			out = append(out, r)
		}
	}
	return out
}
