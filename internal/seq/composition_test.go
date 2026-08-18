package seq

import (
	"math"
	"testing"
)

func TestATContent(t *testing.T) {
	if got := ATContent("ACGT"); math.Abs(got-50) > 0.01 {
		t.Fatalf("ATContent(ACGT) = %f", got)
	}
	if got := ATContent("AATT"); got != 100 {
		t.Fatalf("ATContent(AATT) = %f", got)
	}
}

func TestPurineCount(t *testing.T) {
	if PurineCount("ACGT") != 2 {
		t.Fatalf("purine count = %d", PurineCount("ACGT"))
	}
}

func TestPyrimidineCount(t *testing.T) {
	if PyrimidineCount("ACGT") != 2 {
		t.Fatalf("pyrimidine count = %d", PyrimidineCount("ACGT"))
	}
}

func TestPurinePyrimidineRatio(t *testing.T) {
	ratio := PurinePyrimidineRatio("ACGT")
	if math.Abs(ratio-1.0) > 0.01 {
		t.Fatalf("ratio = %f", ratio)
	}
}

func TestNucleotideFreq(t *testing.T) {
	freq := NucleotideFreq("AACCGG")
	if freq['A'] != 2 || freq['C'] != 2 || freq['G'] != 2 {
		t.Fatalf("freq = %v", freq)
	}
}

func TestDinucleotideFreq(t *testing.T) {
	freq := DinucleotideFreq("ACGT")
	if freq["AC"] != 1 || freq["CG"] != 1 || freq["GT"] != 1 {
		t.Fatalf("freq = %v", freq)
	}
}

func TestLinguisticComplexity(t *testing.T) {
	// "AAAA" is low complexity.
	low := LinguisticComplexity("AAAA", 2)
	// "ACGT" is higher complexity.
	high := LinguisticComplexity("ACGT", 2)
	if high <= low {
		t.Fatalf("ACGT complexity %f should be > AAAA complexity %f", high, low)
	}
}

func TestShannonEntropy(t *testing.T) {
	// Equal composition: entropy should be 2 bits.
	h := ShannonEntropy("ACGT")
	if math.Abs(h-2.0) > 0.01 {
		t.Fatalf("entropy = %f, want 2.0", h)
	}
}

func TestMeltingTemperature(t *testing.T) {
	// Short: ACGT → Tm = 2*(1+1) + 4*(1+1) = 4+8 = 12.
	tm := MeltingTemperature("ACGT")
	if math.Abs(tm-12) > 0.01 {
		t.Fatalf("Tm(ACGT) = %f, want 12", tm)
	}
}
