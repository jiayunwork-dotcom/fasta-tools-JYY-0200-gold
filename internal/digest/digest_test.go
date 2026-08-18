package digest

import "testing"

func TestFindSitesEcoRI(t *testing.T) {
	// EcoRI recognizes GAATTC.
	seq := "AAAGAATTCCCGAATTCAAA"
	enz := CommonEnzymes()["EcoRI"]
	sites, err := FindSites(seq, enz)
	if err != nil {
		t.Fatalf("find sites: %v", err)
	}
	if len(sites) != 2 {
		t.Fatalf("expected 2 sites, got %d", len(sites))
	}
}

func TestDigestEcoRI(t *testing.T) {
	seq := "AAAGAATTCBBBGAATTCCCC"
	enz := CommonEnzymes()["EcoRI"]
	frags, err := Digest(seq, enz)
	if err != nil {
		t.Fatalf("digest: %v", err)
	}
	if len(frags) != 3 {
		t.Fatalf("expected 3 fragments, got %d", len(frags))
	}
	total := 0
	for _, f := range frags {
		total += f.Length
	}
	if total != len(seq) {
		t.Fatalf("total fragment length %d != seq length %d", total, len(seq))
	}
}

func TestDigestNoSite(t *testing.T) {
	seq := "AAAAAAAAA"
	enz := CommonEnzymes()["EcoRI"]
	frags, err := Digest(seq, enz)
	if err != nil {
		t.Fatalf("digest: %v", err)
	}
	if len(frags) != 1 || frags[0].Seq != seq {
		t.Fatalf("expected single fragment of full seq")
	}
}

func TestDigestMulti(t *testing.T) {
	seq := "GAATTCCCCCGGATCC"
	enzymes := []Enzyme{
		CommonEnzymes()["EcoRI"],
		CommonEnzymes()["BamHI"],
	}
	frags, err := DigestMulti(seq, enzymes)
	if err != nil {
		t.Fatalf("digest multi: %v", err)
	}
	if len(frags) < 2 {
		t.Fatalf("expected >= 2 fragments, got %d", len(frags))
	}
}

func TestDigestEmpty(t *testing.T) {
	enz := CommonEnzymes()["EcoRI"]
	_, err := Digest("", enz)
	if err != ErrEmptySeq {
		t.Fatalf("expected ErrEmptySeq, got %v", err)
	}
}

func TestFragmentLengths(t *testing.T) {
	frags := []Fragment{
		{Length: 5},
		{Length: 10},
		{Length: 3},
	}
	lens := FragmentLengths(frags)
	if lens[0] != 10 || lens[1] != 5 || lens[2] != 3 {
		t.Fatalf("not sorted descending: %v", lens)
	}
}
