// Package convert provides format conversion utilities for FASTA records. It
// can emit sequences in various output formats (tabular, JSON lines, multi-line
// FASTA with configurable line width) and filter/transform sequences.
package convert

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"strings"

	"fasta-tools/internal/fasta"
)

// ErrInvalidWidth is returned when line width is <= 0.
var ErrInvalidWidth = errors.New("convert: invalid line width")

// ErrNilWriter is returned when a nil writer is supplied.
var ErrNilWriter = errors.New("convert: nil writer")

// Format enumerates supported output formats.
type Format int

const (
	FormatFASTA  Format = iota // standard multi-line FASTA
	FormatTab                  // tab-separated (header\tsequence)
	FormatJSON                 // JSON Lines (one JSON object per record)
	FormatSingle               // single-line FASTA (no wrapping)
)

// WriteFASTA writes records in standard FASTA format with the given line width.
func WriteFASTA(w io.Writer, records []fasta.Record, lineWidth int) error {
	if w == nil {
		return ErrNilWriter
	}
	if lineWidth <= 0 {
		return ErrInvalidWidth
	}
	for _, r := range records {
		if _, err := fmt.Fprintf(w, ">%s\n", r.Header); err != nil {
			return err
		}
		if err := writeWrapped(w, r.Sequence, lineWidth); err != nil {
			return err
		}
	}
	return nil
}

// writeWrapped writes a sequence string with line breaks every width characters.
func writeWrapped(w io.Writer, seq string, width int) error {
	for i := 0; i < len(seq); i += width {
		end := i + width
		if end > len(seq) {
			end = len(seq)
		}
		if _, err := fmt.Fprintf(w, "%s\n", seq[i:end]); err != nil {
			return err
		}
	}
	if len(seq) == 0 {
		if _, err := fmt.Fprint(w, "\n"); err != nil {
			return err
		}
	}
	return nil
}

// WriteTab writes records in tab-separated format (header\tsequence).
func WriteTab(w io.Writer, records []fasta.Record) error {
	if w == nil {
		return ErrNilWriter
	}
	for _, r := range records {
		if _, err := fmt.Fprintf(w, "%s\t%s\n", r.Header, r.Sequence); err != nil {
			return err
		}
	}
	return nil
}

// jsonRecord is the JSON representation of a FASTA record.
type jsonRecord struct {
	Header   string `json:"header"`
	Sequence string `json:"sequence"`
	Length   int    `json:"length"`
}

// WriteJSON writes records as JSON Lines (one JSON object per line).
func WriteJSON(w io.Writer, records []fasta.Record) error {
	if w == nil {
		return ErrNilWriter
	}
	enc := json.NewEncoder(w)
	enc.SetEscapeHTML(false)
	for _, r := range records {
		jr := jsonRecord{Header: r.Header, Sequence: r.Sequence, Length: len(r.Sequence)}
		if err := enc.Encode(jr); err != nil {
			return err
		}
	}
	return nil
}

// WriteSingle writes records in single-line FASTA (no line wrapping).
func WriteSingle(w io.Writer, records []fasta.Record) error {
	if w == nil {
		return ErrNilWriter
	}
	for _, r := range records {
		if _, err := fmt.Fprintf(w, ">%s\n%s\n", r.Header, r.Sequence); err != nil {
			return err
		}
	}
	return nil
}

// ToUpperCase converts all sequences in records to uppercase.
func ToUpperCase(records []fasta.Record) []fasta.Record {
	out := make([]fasta.Record, len(records))
	for i, r := range records {
		out[i] = fasta.Record{Header: r.Header, Sequence: strings.ToUpper(r.Sequence)}
	}
	return out
}

// ToLowerCase converts all sequences in records to lowercase.
func ToLowerCase(records []fasta.Record) []fasta.Record {
	out := make([]fasta.Record, len(records))
	for i, r := range records {
		out[i] = fasta.Record{Header: r.Header, Sequence: strings.ToLower(r.Sequence)}
	}
	return out
}

// Rename renames records by applying a prefix and sequential numbering.
func Rename(records []fasta.Record, prefix string) []fasta.Record {
	out := make([]fasta.Record, len(records))
	for i, r := range records {
		out[i] = fasta.Record{Header: fmt.Sprintf("%s%d", prefix, i+1), Sequence: r.Sequence}
	}
	return out
}

// FilterByHeader returns records whose header contains the substring.
func FilterByHeader(records []fasta.Record, substr string) []fasta.Record {
	var out []fasta.Record
	for _, r := range records {
		if strings.Contains(r.Header, substr) {
			out = append(out, r)
		}
	}
	return out
}

// RemoveDuplicates removes records with identical sequences, keeping the first.
func RemoveDuplicates(records []fasta.Record) []fasta.Record {
	seen := make(map[string]bool)
	var out []fasta.Record
	for _, r := range records {
		if !seen[r.Sequence] {
			seen[r.Sequence] = true
			out = append(out, r)
		}
	}
	return out
}

// SplitByCount splits records into chunks of at most n records each.
func SplitByCount(records []fasta.Record, n int) [][]fasta.Record {
	if n <= 0 {
		n = 1
	}
	var chunks [][]fasta.Record
	for i := 0; i < len(records); i += n {
		end := i + n
		if end > len(records) {
			end = len(records)
		}
		chunks = append(chunks, records[i:end])
	}
	return chunks
}
