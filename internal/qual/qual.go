// Package qual provides quality scoring utilities for nucleotide sequences.
// It implements Phred quality score conversions, quality filtering, and
// sliding-window quality analysis commonly used in sequencing data QC.
package qual

import (
	"errors"
	"math"
)

// ErrEmptyScores is returned when an empty score slice is provided.
var ErrEmptyScores = errors.New("qual: empty scores")

// ErrInvalidPhred is returned when a Phred score is out of range.
var ErrInvalidPhred = errors.New("qual: phred score out of range [0,93]")

// ErrWindowTooLarge is returned when the window size exceeds the data length.
var ErrWindowTooLarge = errors.New("qual: window larger than data")

// PhredToProb converts a Phred quality score to an error probability.
// Q = -10 * log10(P) → P = 10^(-Q/10)
func PhredToProb(q int) float64 {
	if q < 0 {
		return 1.0
	}
	return math.Pow(10, float64(-q)/10.0)
}

// ProbToPhred converts an error probability to a Phred quality score.
// Returns a clamped value in [0, 93].
func ProbToPhred(p float64) int {
	if p <= 0 {
		return 93
	}
	if p >= 1 {
		return 0
	}
	q := int(-10 * math.Log10(p))
	if q > 93 {
		q = 93
	}
	if q < 0 {
		q = 0
	}
	return q
}

// PhredFromASCII converts a Phred+33 (Sanger/Illumina 1.8+) ASCII character
// to a numeric quality score.
func PhredFromASCII(c byte) int {
	return int(c) - 33
}

// ASCIIFromPhred converts a numeric Phred score to its Phred+33 ASCII encoding.
func ASCIIFromPhred(q int) byte {
	return byte(q + 33)
}

// MeanQuality computes the arithmetic mean of a slice of Phred scores.
func MeanQuality(scores []int) (float64, error) {
	if len(scores) == 0 {
		return 0, ErrEmptyScores
	}
	sum := 0
	for _, s := range scores {
		sum += s
	}
	return float64(sum) / float64(len(scores)), nil
}

// MeanErrorRate computes the mean error rate from Phred scores (average P then
// back to Q). This is more statistically correct than averaging Q values.
func MeanErrorRate(scores []int) (float64, error) {
	if len(scores) == 0 {
		return 0, ErrEmptyScores
	}
	var sumP float64
	for _, q := range scores {
		sumP += PhredToProb(q)
	}
	return sumP / float64(len(scores)), nil
}

// TrimBWA implements the BWA-style quality trimming algorithm. It finds the
// longest suffix of scores where the cumulative sum of (threshold - Q) is
// minimized, effectively trimming low-quality 3' ends.
func TrimBWA(scores []int, threshold int) int {
	if len(scores) == 0 {
		return 0
	}
	n := len(scores)
	bestEnd := n
	maxSum := 0
	runSum := 0
	for i := n - 1; i >= 0; i-- {
		runSum += threshold - scores[i]
		if runSum < 0 {
			runSum = 0
		}
		if runSum > maxSum {
			maxSum = runSum
			bestEnd = i
		}
	}
	return bestEnd
}

// SlidingWindowMean computes the mean quality in a sliding window across
// the scores slice. Returns a slice of means with length (len(scores) - windowSize + 1).
func SlidingWindowMean(scores []int, windowSize int) ([]float64, error) {
	if len(scores) == 0 {
		return nil, ErrEmptyScores
	}
	if windowSize > len(scores) {
		return nil, ErrWindowTooLarge
	}
	if windowSize <= 0 {
		windowSize = 1
	}
	n := len(scores) - windowSize + 1
	means := make([]float64, n)
	// Initial window sum.
	sum := 0
	for i := 0; i < windowSize; i++ {
		sum += scores[i]
	}
	means[0] = float64(sum) / float64(windowSize)
	for i := 1; i < n; i++ {
		sum += scores[i+windowSize-1] - scores[i-1]
		means[i] = float64(sum) / float64(windowSize)
	}
	return means, nil
}

// CountAbove returns the number of scores >= threshold.
func CountAbove(scores []int, threshold int) int {
	n := 0
	for _, s := range scores {
		if s >= threshold {
			n++
		}
	}
	return n
}

// CountBelow returns the number of scores < threshold.
func CountBelow(scores []int, threshold int) int {
	n := 0
	for _, s := range scores {
		if s < threshold {
			n++
		}
	}
	return n
}

// Validate checks that all scores are in the valid Phred range [0, 93].
func Validate(scores []int) error {
	for _, s := range scores {
		if s < 0 || s > 93 {
			return ErrInvalidPhred
		}
	}
	return nil
}

// PercentAbove returns the percentage of positions with quality >= threshold.
func PercentAbove(scores []int, threshold int) float64 {
	if len(scores) == 0 {
		return 0
	}
	return float64(CountAbove(scores, threshold)) / float64(len(scores)) * 100
}
