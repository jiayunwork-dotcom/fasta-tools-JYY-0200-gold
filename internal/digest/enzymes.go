// enzymes.go provides extended enzyme database and helper utilities.
package digest

import "strings"

// EnzymeDB is a searchable collection of restriction enzymes.
type EnzymeDB struct {
	entries map[string]Enzyme
}

// NewEnzymeDB creates a database seeded with the common enzymes.
func NewEnzymeDB() *EnzymeDB {
	db := &EnzymeDB{entries: make(map[string]Enzyme)}
	for name, enz := range CommonEnzymes() {
		db.entries[name] = enz
	}
	return db
}

// Add inserts or updates an enzyme in the database.
func (db *EnzymeDB) Add(e Enzyme) {
	db.entries[e.Name] = e
}

// Get retrieves an enzyme by name (case-insensitive).
func (db *EnzymeDB) Get(name string) (Enzyme, bool) {
	for n, e := range db.entries {
		if strings.EqualFold(n, name) {
			return e, true
		}
	}
	return Enzyme{}, false
}

// SearchBySite returns all enzymes whose recognition site contains the given
// substring.
func (db *EnzymeDB) SearchBySite(sub string) []Enzyme {
	upper := strings.ToUpper(sub)
	var matches []Enzyme
	for _, e := range db.entries {
		if strings.Contains(e.Site, upper) {
			matches = append(matches, e)
		}
	}
	return matches
}

// ListAll returns all enzymes in the database.
func (db *EnzymeDB) ListAll() []Enzyme {
	out := make([]Enzyme, 0, len(db.entries))
	for _, e := range db.entries {
		out = append(out, e)
	}
	return out
}

// Count returns the number of enzymes in the database.
func (db *EnzymeDB) Count() int { return len(db.entries) }

// Remove deletes an enzyme from the database.
func (db *EnzymeDB) Remove(name string) {
	delete(db.entries, name)
}

// FindCutters returns all enzymes that cut the given sequence at least once.
func (db *EnzymeDB) FindCutters(seq string) []Enzyme {
	var cutters []Enzyme
	for _, e := range db.entries {
		sites, err := FindSites(seq, e)
		if err == nil && len(sites) > 0 {
			cutters = append(cutters, e)
		}
	}
	return cutters
}

// FindNonCutters returns all enzymes that do NOT cut the given sequence.
func (db *EnzymeDB) FindNonCutters(seq string) []Enzyme {
	var nonCutters []Enzyme
	for _, e := range db.entries {
		sites, err := FindSites(seq, e)
		if err == nil && len(sites) == 0 {
			nonCutters = append(nonCutters, e)
		}
	}
	return nonCutters
}
