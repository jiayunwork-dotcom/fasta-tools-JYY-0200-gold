// Package orf provides open reading frame (ORF) detection in nucleotide
// sequences. An ORF is defined as a stretch of DNA that begins with a start
// codon (ATG) and ends with a stop codon (TAA, TAG, TGA), with no internal
// stop codons.
//
// The package supports searching all three forward frames and optionally
// the three reverse complement frames.
package orf

import (
	"errors"
	"sort"
	"strings"
)

// ErrEmptySequence is returned when the input sequence is empty.
var ErrEmptySequence = errors.New("orf: empty sequence")

// ORF represents an open reading frame found in a sequence.
type ORF struct {
	Start  int    // 0-based start position in the original sequence
	End    int    // 0-based end position (exclusive, past stop codon)
	Length int    // length in nucleotides
	Frame  int    // reading frame (0, 1, or 2)
	Strand int    // +1 for forward, -1 for reverse complement
	Seq    string // nucleotide sequence of the ORF
}

// AminoAcids returns the length in amino acids (ORF nucleotides / 3).
func (o *ORF) AminoAcids() int { return o.Length / 3 }

// stopCodons contains the standard DNA stop codons.
var stopCodons = map[string]bool{
	"TAA": true, "TAG": true, "TGA": true,
}

// Options configures ORF finding behavior.
type Options struct {
	MinLength    int  // minimum ORF length in nucleotides (default 100)
	BothStrands  bool // search reverse complement as well
	AllowNested  bool // allow ORFs that overlap or are nested
}

// DefaultOptions returns sensible default ORF search options.
func DefaultOptions() Options {
	return Options{MinLength: 100, BothStrands: false, AllowNested: true}
}

// Find locates all ORFs in the given DNA sequence according to the options.
func Find(seq string, opts Options) ([]ORF, error) {
	if len(seq) == 0 {
		return nil, ErrEmptySequence
	}
	upper := strings.ToUpper(seq)
	var orfs []ORF
	// Search forward strand.
	fwd := findInStrand(upper, +1, opts.MinLength)
	orfs = append(orfs, fwd...)
	// Search reverse strand.
	if opts.BothStrands {
		rc := revComp(upper)
		rev := findInStrand(rc, -1, opts.MinLength)
		// Adjust positions to reference the original forward strand.
		for i := range rev {
			origStart := len(seq) - rev[i].End
			origEnd := len(seq) - rev[i].Start
			rev[i].Start = origStart
			rev[i].End = origEnd
		}
		orfs = append(orfs, rev...)
	}
	// Sort by start position.
	sort.Slice(orfs, func(a, b int) bool { return orfs[a].Start < orfs[b].Start })
	return orfs, nil
}

// findInStrand finds ORFs in one strand for all three reading frames.
func findInStrand(seq string, strand int, minLen int) []ORF {
	var orfs []ORF
	for frame := 0; frame < 3; frame++ {
		frameORFs := findInFrame(seq, frame, strand, minLen)
		orfs = append(orfs, frameORFs...)
	}
	return orfs
}

// findInFrame finds all ORFs in a specific frame of the sequence.
func findInFrame(seq string, frame, strand, minLen int) []ORF {
	var orfs []ORF
	n := len(seq)
	i := frame
	for i+3 <= n {
		codon := seq[i : i+3]
		if codon != "ATG" {
			i += 3
			continue
		}
		// Found start codon; scan for stop.
		start := i
		j := i + 3
		for j+3 <= n {
			c := seq[j : j+3]
			if stopCodons[c] {
				orfLen := j + 3 - start
				if orfLen >= minLen {
					orfs = append(orfs, ORF{
						Start:  start,
						End:    j + 3,
						Length: orfLen,
						Frame:  frame,
						Strand: strand,
						Seq:    seq[start : j+3],
					})
				}
				i = j + 3
				goto next
			}
			j += 3
		}
		// No stop codon found; move past this start.
		i += 3
	next:
	}
	return orfs
}

// Longest returns the longest ORF from the list, or nil if empty.
func Longest(orfs []ORF) *ORF {
	if len(orfs) == 0 {
		return nil
	}
	best := &orfs[0]
	for i := 1; i < len(orfs); i++ {
		if orfs[i].Length > best.Length {
			best = &orfs[i]
		}
	}
	return best
}

// FilterByLength returns only ORFs at least minLen nucleotides long.
func FilterByLength(orfs []ORF, minLen int) []ORF {
	var out []ORF
	for _, o := range orfs {
		if o.Length >= minLen {
			out = append(out, o)
		}
	}
	return out
}

// CountByFrame returns the number of ORFs in each reading frame.
func CountByFrame(orfs []ORF) [3]int {
	var counts [3]int
	for _, o := range orfs {
		if o.Frame >= 0 && o.Frame < 3 {
			counts[o.Frame]++
		}
	}
	return counts
}

// revComp returns the reverse complement of a DNA sequence.
func revComp(s string) string {
	out := make([]byte, len(s))
	for i := range s {
		var c byte
		switch s[len(s)-1-i] {
		case 'A':
			c = 'T'
		case 'T':
			c = 'A'
		case 'C':
			c = 'G'
		case 'G':
			c = 'C'
		default:
			c = 'N'
		}
		out[i] = c
	}
	return string(out)
}
