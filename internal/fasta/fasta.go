package fasta

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

// Record holds a single FASTA entry.
type Record struct {
	Header   string
	Sequence string
}

var validSeq = func() map[byte]bool {
	m := make(map[byte]bool)
	for _, c := range "ACGTUacgtuN" {
		m[byte(c)] = true
	}
	return m
}()

// Parse reads FASTA records from r.
// Every header line must start with '>'. Sequence characters are
// restricted to ACGTUacgtuN. A bad character or sequence data before any
// header causes an error.
func Parse(r io.Reader) ([]Record, error) {
	var records []Record

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)

	var cur *Record
	var seqBuilder strings.Builder

	flush := func() error {
		if cur == nil {
			return nil
		}
		seq := seqBuilder.String()
		for i := 0; i < len(seq); i++ {
			if !validSeq[seq[i]] {
				return fmt.Errorf("invalid character %q in sequence of %q", seq[i], cur.Header)
			}
		}
		cur.Sequence = seq
		records = append(records, *cur)
		return nil
	}

	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := scanner.Text()
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			if err := flush(); err != nil {
				return nil, err
			}
			cur = &Record{Header: strings.TrimSpace(strings.TrimPrefix(line, ">"))}
			seqBuilder.Reset()
			continue
		}
		if cur == nil {
			return nil, fmt.Errorf("sequence data before any header at line %d", lineNo)
		}
		seqBuilder.WriteString(strings.TrimSpace(line))
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	if err := flush(); err != nil {
		return nil, err
	}
	if len(records) == 0 {
		return nil, fmt.Errorf("no FASTA records found")
	}
	return records, nil
}
