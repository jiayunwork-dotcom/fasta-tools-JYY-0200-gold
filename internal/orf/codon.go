// codon.go extends the orf package with codon analysis utilities that operate
// on identified ORFs: codon bias, start/stop codon statistics, and reading
// frame analysis.
package orf

import (
	"strings"
)

// CodonBias computes the codon usage frequencies within an ORF. Returns a map
// from codon (3-letter DNA string) to count.
func CodonBias(o ORF) map[string]int {
	upper := strings.ToUpper(o.Seq)
	counts := make(map[string]int)
	for i := 0; i+3 <= len(upper); i += 3 {
		codon := upper[i : i+3]
		counts[codon]++
	}
	return counts
}

// StartCodonContext extracts the nucleotide context around the start codon of
// an ORF. Returns up to `flanking` bases on each side of the ATG.
func StartCodonContext(seq string, o ORF, flanking int) string {
	start := o.Start - flanking
	if start < 0 {
		start = 0
	}
	end := o.Start + 3 + flanking
	if end > len(seq) {
		end = len(seq)
	}
	return seq[start:end]
}

// StopCodonType returns the stop codon at the end of the ORF (TAA, TAG, TGA).
func StopCodonType(o ORF) string {
	if len(o.Seq) < 3 {
		return ""
	}
	return strings.ToUpper(o.Seq[len(o.Seq)-3:])
}

// GCAtPosition computes the GC percentage at each codon position (1st, 2nd, 3rd)
// within the ORF.
func GCAtPosition(o ORF) [3]float64 {
	upper := strings.ToUpper(o.Seq)
	var counts [3]int
	var gc [3]int
	for i := 0; i+3 <= len(upper); i += 3 {
		for p := 0; p < 3; p++ {
			counts[p]++
			if upper[i+p] == 'G' || upper[i+p] == 'C' {
				gc[p]++
			}
		}
	}
	var result [3]float64
	for p := 0; p < 3; p++ {
		if counts[p] > 0 {
			result[p] = float64(gc[p]) / float64(counts[p]) * 100
		}
	}
	return result
}

// FrameDistribution calculates what fraction of the total ORF nucleotides
// falls in each reading frame (0, 1, 2).
func FrameDistribution(orfs []ORF) [3]float64 {
	var totals [3]int
	sum := 0
	for _, o := range orfs {
		totals[o.Frame] += o.Length
		sum += o.Length
	}
	var dist [3]float64
	if sum > 0 {
		for i := 0; i < 3; i++ {
			dist[i] = float64(totals[i]) / float64(sum) * 100
		}
	}
	return dist
}

// AvgORFLength computes the average ORF length in nucleotides.
func AvgORFLength(orfs []ORF) float64 {
	if len(orfs) == 0 {
		return 0
	}
	total := 0
	for _, o := range orfs {
		total += o.Length
	}
	return float64(total) / float64(len(orfs))
}

// CodingDensity computes the fraction of the sequence covered by ORFs.
func CodingDensity(seqLen int, orfs []ORF) float64 {
	if seqLen == 0 {
		return 0
	}
	covered := 0
	for _, o := range orfs {
		covered += o.Length
	}
	ratio := float64(covered) / float64(seqLen)
	if ratio > 1.0 {
		ratio = 1.0
	}
	return ratio
}
