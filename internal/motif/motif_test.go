package motif

import (
	"math"
	"testing"
)

func TestFindExact(t *testing.T) {
	matches, err := FindExact("ACGTACGTACGT", "ACGT")
	if err != nil {
		t.Fatalf("find: %v", err)
	}
	if len(matches) != 3 {
		t.Fatalf("expected 3 matches, got %d", len(matches))
	}
	if matches[0].Start != 0 || matches[1].Start != 4 || matches[2].Start != 8 {
		t.Fatalf("wrong positions: %v", matches)
	}
}

func TestFindExactCaseInsensitive(t *testing.T) {
	matches, err := FindExact("acgtACGT", "ACGT")
	if err != nil {
		t.Fatalf("find: %v", err)
	}
	if len(matches) != 2 {
		t.Fatalf("expected 2 matches, got %d", len(matches))
	}
}

func TestFindExactEmpty(t *testing.T) {
	_, err := FindExact("ACGT", "")
	if err != ErrEmptyPattern {
		t.Fatalf("expected ErrEmptyPattern, got %v", err)
	}
}

func TestFindIUPAC(t *testing.T) {
	// R matches A or G.
	matches, err := FindIUPAC("ACGT", "RCG")
	if err != nil {
		t.Fatalf("find: %v", err)
	}
	if len(matches) != 1 || matches[0].Start != 0 {
		t.Fatalf("expected 1 match at 0, got %v", matches)
	}
}

func TestFindIUPACMultiple(t *testing.T) {
	// N matches any base.
	matches, err := FindIUPAC("ACGT", "N")
	if err != nil {
		t.Fatalf("find: %v", err)
	}
	if len(matches) != 4 {
		t.Fatalf("expected 4 matches, got %d", len(matches))
	}
}

func TestFindIUPACInvalid(t *testing.T) {
	_, err := FindIUPAC("ACGT", "Z")
	if err != ErrInvalidIUPAC {
		t.Fatalf("expected ErrInvalidIUPAC, got %v", err)
	}
}

func TestPWMScan(t *testing.T) {
	// Create a simple PWM that strongly prefers "ACGT".
	freq := [][4]float64{
		{10, 0, 0, 0}, // A
		{0, 10, 0, 0}, // C
		{0, 0, 10, 0}, // G
		{0, 0, 0, 10}, // T
	}
	pwm := NewPWM("test", freq, 0.1)
	matches := pwm.Scan("NNACGTNNN", 0)
	if len(matches) == 0 {
		t.Fatal("expected at least one PWM match")
	}
	// The best match should be at position 2 (ACGT).
	best := matches[0]
	for _, m := range matches[1:] {
		if m.Score > best.Score {
			best = m
		}
	}
	if best.Start != 2 {
		t.Fatalf("best match at %d, want 2", best.Start)
	}
}

func TestPWMMaxScore(t *testing.T) {
	freq := [][4]float64{
		{5, 5, 5, 5},
	}
	pwm := NewPWM("uniform", freq, 0.1)
	// Uniform frequencies should give ~0 log-odds against uniform background.
	if math.Abs(pwm.MaxScore()) > 0.5 {
		t.Fatalf("max score should be ~0, got %f", pwm.MaxScore())
	}
}

func TestCountOccurrences(t *testing.T) {
	n := CountOccurrences("AAAA", "AA")
	if n != 3 {
		t.Fatalf("count = %d, want 3", n)
	}
}
