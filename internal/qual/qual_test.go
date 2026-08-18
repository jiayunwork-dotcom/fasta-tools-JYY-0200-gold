package qual

import (
	"math"
	"testing"
)

func TestPhredToProb(t *testing.T) {
	cases := []struct {
		q    int
		want float64
	}{
		{0, 1.0},
		{10, 0.1},
		{20, 0.01},
		{30, 0.001},
	}
	for _, c := range cases {
		got := PhredToProb(c.q)
		if math.Abs(got-c.want) > 1e-9 {
			t.Fatalf("PhredToProb(%d) = %f, want %f", c.q, got, c.want)
		}
	}
}

func TestProbToPhred(t *testing.T) {
	if ProbToPhred(0.01) != 20 {
		t.Fatalf("ProbToPhred(0.01) = %d", ProbToPhred(0.01))
	}
	if ProbToPhred(0) != 93 {
		t.Fatalf("ProbToPhred(0) = %d", ProbToPhred(0))
	}
	if ProbToPhred(1.0) != 0 {
		t.Fatalf("ProbToPhred(1.0) = %d", ProbToPhred(1.0))
	}
}

func TestPhredASCII(t *testing.T) {
	for q := 0; q < 42; q++ {
		c := ASCIIFromPhred(q)
		got := PhredFromASCII(c)
		if got != q {
			t.Fatalf("round-trip %d -> %c -> %d", q, c, got)
		}
	}
}

func TestMeanQuality(t *testing.T) {
	scores := []int{30, 30, 30, 30}
	m, err := MeanQuality(scores)
	if err != nil {
		t.Fatalf("mean: %v", err)
	}
	if m != 30 {
		t.Fatalf("mean = %f, want 30", m)
	}
}

func TestMeanQualityEmpty(t *testing.T) {
	_, err := MeanQuality(nil)
	if err != ErrEmptyScores {
		t.Fatalf("expected ErrEmptyScores, got %v", err)
	}
}

func TestTrimBWA(t *testing.T) {
	// High quality then low quality tail.
	scores := []int{30, 30, 30, 30, 5, 5, 5}
	trimPos := TrimBWA(scores, 15)
	if trimPos > 4 {
		t.Fatalf("trim pos = %d, expected <= 4", trimPos)
	}
}

func TestSlidingWindowMean(t *testing.T) {
	scores := []int{10, 20, 30, 40, 50}
	means, err := SlidingWindowMean(scores, 3)
	if err != nil {
		t.Fatalf("sliding: %v", err)
	}
	// First window: (10+20+30)/3 = 20.
	if math.Abs(means[0]-20) > 0.01 {
		t.Fatalf("means[0] = %f, want 20", means[0])
	}
	if len(means) != 3 {
		t.Fatalf("len = %d, want 3", len(means))
	}
}

func TestSlidingWindowTooLarge(t *testing.T) {
	_, err := SlidingWindowMean([]int{1, 2}, 5)
	if err != ErrWindowTooLarge {
		t.Fatalf("expected ErrWindowTooLarge, got %v", err)
	}
}

func TestCountAboveBelow(t *testing.T) {
	scores := []int{10, 20, 30, 40, 50}
	if CountAbove(scores, 30) != 3 {
		t.Fatalf("above 30 = %d", CountAbove(scores, 30))
	}
	if CountBelow(scores, 30) != 2 {
		t.Fatalf("below 30 = %d", CountBelow(scores, 30))
	}
}

func TestValidate(t *testing.T) {
	if err := Validate([]int{0, 40, 93}); err != nil {
		t.Fatalf("valid scores: %v", err)
	}
	if err := Validate([]int{-1}); err != ErrInvalidPhred {
		t.Fatalf("expected ErrInvalidPhred, got %v", err)
	}
	if err := Validate([]int{94}); err != ErrInvalidPhred {
		t.Fatalf("expected ErrInvalidPhred, got %v", err)
	}
}

func TestPercentAbove(t *testing.T) {
	scores := []int{10, 20, 30, 40}
	pct := PercentAbove(scores, 25)
	if math.Abs(pct-50) > 0.01 {
		t.Fatalf("percent above = %f, want 50", pct)
	}
}
