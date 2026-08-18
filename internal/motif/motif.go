// Package motif provides motif and pattern searching in nucleotide sequences.
// It supports exact substring matching, IUPAC ambiguity code matching, and
// position weight matrix (PWM) scanning for identifying binding sites or
// regulatory elements.
package motif

import (
	"errors"
	"math"
	"strings"
)

// ErrEmptyPattern is returned when an empty pattern is provided.
var ErrEmptyPattern = errors.New("motif: empty pattern")

// ErrInvalidIUPAC is returned when an unknown IUPAC code is encountered.
var ErrInvalidIUPAC = errors.New("motif: invalid IUPAC code")

// Match represents a motif occurrence in a sequence.
type Match struct {
	Start   int     // 0-based start position
	End     int     // 0-based end position (exclusive)
	Seq     string  // matched subsequence
	Score   float64 // match score (1.0 for exact, variable for PWM)
	Strand  int     // +1 for forward, -1 for reverse complement
}

// iupac maps IUPAC ambiguity codes to the set of matching bases.
var iupac = map[byte]string{
	'A': "A", 'C': "C", 'G': "G", 'T': "T", 'U': "U",
	'R': "AG", 'Y': "CT", 'S': "GC", 'W': "AT", 'K': "GT", 'M': "AC",
	'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG",
	'N': "ACGT",
}

// FindExact finds all exact occurrences of pattern in seq (case insensitive).
func FindExact(seq, pattern string) ([]Match, error) {
	if len(pattern) == 0 {
		return nil, ErrEmptyPattern
	}
	upper := strings.ToUpper(seq)
	pat := strings.ToUpper(pattern)
	var matches []Match
	for i := 0; i+len(pat) <= len(upper); i++ {
		if upper[i:i+len(pat)] == pat {
			matches = append(matches, Match{
				Start:  i,
				End:    i + len(pat),
				Seq:    seq[i : i+len(pat)],
				Score:  1.0,
				Strand: 1,
			})
		}
	}
	return matches, nil
}

// FindIUPAC finds all positions where the sequence matches an IUPAC pattern.
// The pattern uses standard IUPAC ambiguity codes (case insensitive).
func FindIUPAC(seq, pattern string) ([]Match, error) {
	if len(pattern) == 0 {
		return nil, ErrEmptyPattern
	}
	upper := strings.ToUpper(seq)
	pat := strings.ToUpper(pattern)
	// Validate pattern.
	for i := range pat {
		if _, ok := iupac[pat[i]]; !ok {
			return nil, ErrInvalidIUPAC
		}
	}
	var matches []Match
	for i := 0; i+len(pat) <= len(upper); i++ {
		if matchIUPAC(upper[i:i+len(pat)], pat) {
			matches = append(matches, Match{
				Start:  i,
				End:    i + len(pat),
				Seq:    seq[i : i+len(pat)],
				Score:  1.0,
				Strand: 1,
			})
		}
	}
	return matches, nil
}

// matchIUPAC checks if sub matches the IUPAC pattern.
func matchIUPAC(sub, pat string) bool {
	for i := range pat {
		allowed := iupac[pat[i]]
		if !strings.ContainsRune(allowed, rune(sub[i])) {
			return false
		}
	}
	return true
}

// PWM is a position weight matrix for motif scoring.
// Rows are positions (length of the motif); columns are bases (A, C, G, T).
type PWM struct {
	Matrix []PWMRow // one row per position
	Name   string
}

// PWMRow holds the log-odds scores for one position.
type PWMRow struct {
	A, C, G, T float64
}

// NewPWM creates a PWM from a frequency matrix. Each row gives counts for
// A, C, G, T. A pseudocount is added and log-odds are computed against a
// uniform background (0.25 each).
func NewPWM(name string, freq [][4]float64, pseudo float64) *PWM {
	rows := make([]PWMRow, len(freq))
	for i, f := range freq {
		total := f[0] + f[1] + f[2] + f[3] + 4*pseudo
		rows[i] = PWMRow{
			A: math.Log2((f[0] + pseudo) / total / 0.25),
			C: math.Log2((f[1] + pseudo) / total / 0.25),
			G: math.Log2((f[2] + pseudo) / total / 0.25),
			T: math.Log2((f[3] + pseudo) / total / 0.25),
		}
	}
	return &PWM{Matrix: rows, Name: name}
}

// Len returns the length (number of positions) of the PWM.
func (p *PWM) Len() int { return len(p.Matrix) }

// ScoreAt scores the subsequence at position i in seq against the PWM.
func (p *PWM) ScoreAt(seq string, i int) float64 {
	score := 0.0
	for k, row := range p.Matrix {
		pos := i + k
		if pos >= len(seq) {
			return math.Inf(-1)
		}
		switch seq[pos] {
		case 'A', 'a':
			score += row.A
		case 'C', 'c':
			score += row.C
		case 'G', 'g':
			score += row.G
		case 'T', 't', 'U', 'u':
			score += row.T
		default:
			score += math.Min(math.Min(row.A, row.C), math.Min(row.G, row.T))
		}
	}
	return score
}

// MaxScore returns the maximum possible score of this PWM.
func (p *PWM) MaxScore() float64 {
	s := 0.0
	for _, row := range p.Matrix {
		s += math.Max(math.Max(row.A, row.C), math.Max(row.G, row.T))
	}
	return s
}

// Scan finds all positions in seq that score above the threshold.
func (p *PWM) Scan(seq string, threshold float64) []Match {
	upper := strings.ToUpper(seq)
	var matches []Match
	for i := 0; i+p.Len() <= len(upper); i++ {
		sc := p.ScoreAt(upper, i)
		if sc >= threshold {
			matches = append(matches, Match{
				Start:  i,
				End:    i + p.Len(),
				Seq:    seq[i : i+p.Len()],
				Score:  sc,
				Strand: 1,
			})
		}
	}
	return matches
}

// CountOccurrences counts the number of times pattern appears in seq (exact).
func CountOccurrences(seq, pattern string) int {
	matches, err := FindExact(seq, pattern)
	if err != nil {
		return 0
	}
	return len(matches)
}
