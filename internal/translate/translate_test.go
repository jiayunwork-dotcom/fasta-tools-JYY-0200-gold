package translate

import "testing"

func TestDNAToRNA(t *testing.T) {
	if got := DNAToRNA("ACGT"); got != "ACGU" {
		t.Fatalf("DNAToRNA(ACGT) = %q", got)
	}
	if got := DNAToRNA("acgt"); got != "acgu" {
		t.Fatalf("DNAToRNA(acgt) = %q", got)
	}
}

func TestRNAToDNA(t *testing.T) {
	if got := RNAToDNA("ACGU"); got != "ACGT" {
		t.Fatalf("RNAToDNA(ACGU) = %q", got)
	}
}

func TestTranslateRNA(t *testing.T) {
	// AUG = M, UUU = F, UAA = *
	prot, err := TranslateRNA("AUGUUUUAA")
	if err != nil {
		t.Fatalf("translate: %v", err)
	}
	if prot != "MF*" {
		t.Fatalf("prot = %q, want MF*", prot)
	}
}

func TestTranslateDNA(t *testing.T) {
	// ATG -> AUG = M
	prot, err := TranslateDNA("ATGTTTCCC")
	if err != nil {
		t.Fatalf("translate: %v", err)
	}
	if prot != "MFP" {
		t.Fatalf("prot = %q, want MFP", prot)
	}
}

func TestTranslateEmpty(t *testing.T) {
	_, err := TranslateRNA("")
	if err != ErrEmptyInput {
		t.Fatalf("expected ErrEmptyInput, got %v", err)
	}
}

func TestTranslateFrame(t *testing.T) {
	// Frame 0: ATG TTT = MF
	// Frame 1: TGT TTC = CF (from ATGTTTCCC -> frame1 = TGTTTCCC -> TGU UUC = CF)
	dna := "ATGTTTCCC"
	p0, _ := TranslateFrame(dna, 0)
	if p0 != "MFP" {
		t.Fatalf("frame0 = %q, want MFP", p0)
	}
	p1, _ := TranslateFrame(dna, 1)
	if len(p1) == 0 {
		t.Fatal("frame1 should produce result")
	}
}

func TestThreeFrame(t *testing.T) {
	frames, err := ThreeFrame("ATGCCCGGGAAA")
	if err != nil {
		t.Fatalf("three frame: %v", err)
	}
	for i, f := range frames {
		if len(f) == 0 {
			t.Fatalf("frame %d is empty", i)
		}
	}
}

func TestSixFrame(t *testing.T) {
	frames, err := SixFrame("ATGCCCGGGAAATTT")
	if err != nil {
		t.Fatalf("six frame: %v", err)
	}
	for i, f := range frames {
		if len(f) == 0 {
			t.Fatalf("frame %d is empty", i)
		}
	}
}

func TestCodonUsage(t *testing.T) {
	counts, err := CodonUsage("ATGATGATG")
	if err != nil {
		t.Fatalf("codon usage: %v", err)
	}
	if counts["AUG"] != 3 {
		t.Fatalf("AUG count = %d, want 3", counts["AUG"])
	}
}

func TestLookupCodon(t *testing.T) {
	aa, ok := LookupCodon("AUG")
	if !ok || aa != 'M' {
		t.Fatalf("AUG -> %c, %v", aa, ok)
	}
	_, ok = LookupCodon("XYZ")
	if ok {
		t.Fatal("XYZ should not be found")
	}
}
