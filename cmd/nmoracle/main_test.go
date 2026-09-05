package main

import (
	"testing"

	"github.com/rainoffallingstar/bamdriver/pkg/bamnative"
)

func TestScoreAlignmentAppliesStrandAwareBisulfiteConversions(t *testing.T) {
	testCases := []struct {
		name      string
		flags     uint16
		reference string
		sequence  string
		operation byte
		wantNM    int64
	}{
		{
			name:      "forward C-to-T conversion is ignored",
			reference: "C",
			sequence:  "T",
			operation: bamnative.CigarMatch,
			wantNM:    0,
		},
		{
			name:      "forward G-to-A remains a mismatch",
			reference: "G",
			sequence:  "A",
			operation: bamnative.CigarMatch,
			wantNM:    1,
		},
		{
			name:      "reverse G-to-A conversion is ignored",
			flags:     bamnative.FlagReverse,
			reference: "G",
			sequence:  "A",
			operation: bamnative.CigarMismatch,
			wantNM:    0,
		},
		{
			name:      "reverse C-to-T remains a mismatch",
			flags:     bamnative.FlagReverse,
			reference: "C",
			sequence:  "T",
			operation: bamnative.CigarMismatch,
			wantNM:    1,
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			record := &bamnative.Record{
				RefID: 0,
				Flags: testCase.flags,
				Seq:   testCase.sequence,
				Cigar: []bamnative.CigarOp{{Op: testCase.operation, Len: 1}},
			}
			score, err := scoreAlignment(record, []byte(testCase.reference))
			if err != nil {
				t.Fatalf("scoreAlignment: %v", err)
			}
			if score.XenofilxBisulfiteNM != testCase.wantNM {
				t.Fatalf("bisulfite NM = %d, want %d", score.XenofilxBisulfiteNM, testCase.wantNM)
			}
		})
	}
}

func TestScoreAlignmentUsesBismarkGenomeConversionContext(t *testing.T) {
	testCases := []struct {
		name      string
		flags     uint16
		context   string
		reference string
		sequence  string
		wantNM    int64
	}{
		{
			name:      "GA context overrides forward flag",
			context:   "GA",
			reference: "G",
			sequence:  "A",
			wantNM:    0,
		},
		{
			name:      "CT context overrides reverse flag",
			flags:     bamnative.FlagReverse,
			context:   "CT",
			reference: "C",
			sequence:  "T",
			wantNM:    0,
		},
		{
			name:      "invalid context falls back to flag",
			context:   "invalid",
			reference: "C",
			sequence:  "T",
			wantNM:    0,
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			record := &bamnative.Record{
				RefID: 0,
				Flags: testCase.flags,
				Seq:   testCase.sequence,
				Cigar: []bamnative.CigarOp{{Op: bamnative.CigarMatch, Len: 1}},
				Aux: []*bamnative.AuxField{
					{Tag: "XG", Type: bamnative.AuxTypeString, Value: testCase.context},
				},
			}
			score, err := scoreAlignment(record, []byte(testCase.reference))
			if err != nil {
				t.Fatalf("scoreAlignment: %v", err)
			}
			if score.XenofilxBisulfiteNM != testCase.wantNM {
				t.Fatalf("bisulfite NM = %d, want %d", score.XenofilxBisulfiteNM, testCase.wantNM)
			}
		})
	}
}
