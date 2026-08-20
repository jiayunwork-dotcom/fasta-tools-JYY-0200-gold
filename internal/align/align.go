// Package align implements pairwise sequence alignment algorithms. It provides
// both global alignment (Needleman-Wunsch) and local alignment (Smith-Waterman)
// using dynamic programming with configurable scoring matrices.
package align

import "errors"

// ErrEmptySequence is returned when one or both input sequences are empty.
var ErrEmptySequence = errors.New("align: empty sequence")

// ErrInvalidGap is returned when gap penalties are positive (must be ≤ 0).
var ErrInvalidGap = errors.New("align: gap penalty must be <= 0")

// Scoring holds the alignment scoring parameters.
type Scoring struct {
	Match    int // reward for matching bases (positive)
	Mismatch int // penalty for mismatches (negative)
	GapOpen  int // penalty for opening a gap (negative)
	GapExt   int // penalty for extending a gap (negative)
}

// DefaultScoring returns a commonly used scoring matrix.
func DefaultScoring() Scoring {
	return Scoring{Match: 2, Mismatch: -1, GapOpen: -2, GapExt: -1}
}

// Result holds the outcome of an alignment.
type Result struct {
	AlignedA string // aligned first sequence (with gaps as '-')
	AlignedB string // aligned second sequence (with gaps as '-')
	Score    int    // alignment score
	Identity float64 // fraction of positions that match [0,1]
	Length   int    // length of the alignment
}

// max3 returns the maximum of three integers and its index (0, 1, or 2).
func max3(a, b, c int) (int, int) {
	m, idx := a, 0
	if b > m {
		m, idx = b, 1
	}
	if c > m {
		m, idx = c, 2
	}
	return m, idx
}

// max4 returns the maximum of four integers.
func max4(a, b, c, d int) int {
	m := a
	if b > m {
		m = b
	}
	if c > m {
		m = c
	}
	if d > m {
		m = d
	}
	return m
}

// Global performs Needleman-Wunsch global alignment between sequences a and b.
// It uses simple linear gap penalty (GapOpen used for each gap position).
func Global(a, b string, sc Scoring) (*Result, error) {
	if len(a) == 0 || len(b) == 0 {
		return nil, ErrEmptySequence
	}
	if sc.GapOpen > 0 || sc.GapExt > 0 {
		return nil, ErrInvalidGap
	}
	m, n := len(a), len(b)
	gap := sc.GapOpen

	// DP matrix.
	dp := make([][]int, m+1)
	for i := range dp {
		dp[i] = make([]int, n+1)
	}
	for i := 1; i <= m; i++ {
		dp[i][0] = dp[i-1][0] + gap
	}
	for j := 1; j <= n; j++ {
		dp[0][j] = dp[0][j-1] + gap
	}
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			s := sc.Mismatch
			if a[i-1] == b[j-1] {
				s = sc.Match
			}
			diag := dp[i-1][j-1] + s
			up := dp[i-1][j] + gap
			left := dp[i][j-1] + gap
			v, _ := max3(diag, up, left)
			dp[i][j] = v
		}
	}

	// Traceback.
	var alignA, alignB []byte
	i, j := m, n
	for i > 0 || j > 0 {
		if i > 0 && j > 0 {
			s := sc.Mismatch
			if a[i-1] == b[j-1] {
				s = sc.Match
			}
			if dp[i][j] == dp[i-1][j-1]+s {
				alignA = append(alignA, a[i-1])
				alignB = append(alignB, b[j-1])
				i--
				j--
				continue
			}
		}
		if i > 0 && dp[i][j] == dp[i-1][j]+gap {
			alignA = append(alignA, a[i-1])
			alignB = append(alignB, '-')
			i--
		} else {
			alignA = append(alignA, '-')
			alignB = append(alignB, b[j-1])
			j--
		}
	}
	// Reverse.
	reverse(alignA)
	reverse(alignB)

	res := &Result{
		AlignedA: string(alignA),
		AlignedB: string(alignB),
		Score:    dp[m][n],
		Length:   len(alignA),
	}
	res.Identity = computeIdentity(alignA, alignB)
	return res, nil
}

// Local performs Smith-Waterman local alignment between sequences a and b.
func Local(a, b string, sc Scoring) (*Result, error) {
	if len(a) == 0 || len(b) == 0 {
		return nil, ErrEmptySequence
	}
	if sc.GapOpen > 0 || sc.GapExt > 0 {
		return nil, ErrInvalidGap
	}
	m, n := len(a), len(b)
	gap := sc.GapOpen

	dp := make([][]int, m+1)
	for i := range dp {
		dp[i] = make([]int, n+1)
	}
	bestScore, bestI, bestJ := 0, 0, 0
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			s := sc.Mismatch
			if a[i-1] == b[j-1] {
				s = sc.Match
			}
			v := max4(0, dp[i-1][j-1]+s, dp[i-1][j]+gap, dp[i][j-1]+gap)
			dp[i][j] = v
			if v > bestScore {
				bestScore = v
				bestI = i
				bestJ = j
			}
		}
	}
	if bestScore == 0 {
		return &Result{Score: 0}, nil
	}

	// Traceback from best cell until we hit 0.
	var alignA, alignB []byte
	i, j := bestI, bestJ
	for i > 0 && j > 0 && dp[i][j] > 0 {
		s := sc.Mismatch
		if a[i-1] == b[j-1] {
			s = sc.Match
		}
		if dp[i][j] == dp[i-1][j-1]+s {
			alignA = append(alignA, a[i-1])
			alignB = append(alignB, b[j-1])
			i--
			j--
		} else if dp[i][j] == dp[i-1][j]+gap {
			alignA = append(alignA, a[i-1])
			alignB = append(alignB, '-')
			i--
		} else {
			alignA = append(alignA, '-')
			alignB = append(alignB, b[j-1])
			j--
		}
	}
	reverse(alignA)
	reverse(alignB)

	res := &Result{
		AlignedA: string(alignA),
		AlignedB: string(alignB),
		Score:    bestScore,
		Length:   len(alignA),
	}
	res.Identity = computeIdentity(alignA, alignB)
	return res, nil
}

// reverse a byte slice in place.
func reverse(b []byte) {
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
}

// computeIdentity calculates the fraction of matching positions.
func computeIdentity(a, b []byte) float64 {
	if len(a) == 0 {
		return 0
	}
	match := 0
	for i := range a {
		if a[i] == b[i] && a[i] != '-' {
			match++
		}
	}
	return float64(match) / float64(len(a))
}

// HammingDistance returns the number of positions where a and b differ.
// Both sequences must have the same length.
func HammingDistance(a, b string) (int, error) {
	if len(a) != len(b) {
		return 0, errors.New("align: sequences must have equal length for Hamming distance")
	}
	d := 0
	for i := range a {
		if a[i] != b[i] {
			d++
		}
	}
	return d, nil
}

// EditDistance computes the Levenshtein edit distance between two sequences.
func EditDistance(a, b string) int {
	m, n := len(a), len(b)
	if m == 0 {
		return n
	}
	if n == 0 {
		return m
	}
	prev := make([]int, n+1)
	curr := make([]int, n+1)
	for j := 0; j <= n; j++ {
		prev[j] = j
	}
	for i := 1; i <= m; i++ {
		curr[0] = i
		for j := 1; j <= n; j++ {
			cost := 1
			if a[i-1] == b[j-1] {
				cost = 0
			}
			ins := curr[j-1] + 1
			del := prev[j] + 1
			sub := prev[j-1] + cost
			mn := ins
			if del < mn {
				mn = del
			}
			if sub < mn {
				mn = sub
			}
			curr[j] = mn
		}
		prev, curr = curr, prev
	}
	return prev[n]
}
