package main

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"

	"github.com/rainoffallingstar/bamdriver/pkg/bamnative"
)

func TestRoundTripAndComparePreserveCanonicalRecords(t *testing.T) {
	temporaryDirectory := t.TempDir()
	inputPath := filepath.Join(temporaryDirectory, "input.bam")
	outputPath := filepath.Join(temporaryDirectory, "output.bam")
	roundTripReportPath := filepath.Join(temporaryDirectory, "roundtrip.json")
	comparisonReportPath := filepath.Join(temporaryDirectory, "comparison.json")

	writeFixtureBAM(t, inputPath)
	if err := runRoundTrip([]string{
		"--input", inputPath,
		"--output", outputPath,
		"--report", roundTripReportPath,
	}); err != nil {
		t.Fatalf("runRoundTrip: %v", err)
	}
	if err := runComparison([]string{
		"--left", inputPath,
		"--right", outputPath,
		"--report", comparisonReportPath,
	}); err != nil {
		t.Fatalf("runComparison: %v", err)
	}

	var report comparisonReport
	readReport(t, comparisonReportPath, &report)
	if !report.Equal {
		t.Fatalf("comparison report = %+v, want equality", report)
	}
	if report.RecordsCompared != 2 {
		t.Fatalf("records compared = %d, want 2", report.RecordsCompared)
	}
}

func TestComparisonReportsRecordMutation(t *testing.T) {
	temporaryDirectory := t.TempDir()
	leftPath := filepath.Join(temporaryDirectory, "left.bam")
	rightPath := filepath.Join(temporaryDirectory, "right.bam")
	reportPath := filepath.Join(temporaryDirectory, "comparison.json")

	writeFixtureBAM(t, leftPath)
	writeMutatedFixtureBAM(t, rightPath)
	comparisonError := runComparison([]string{
		"--left", leftPath,
		"--right", rightPath,
		"--report", reportPath,
	})
	if comparisonError == nil {
		t.Fatal("runComparison succeeded for a mutated BAM")
	}

	var report comparisonReport
	readReport(t, reportPath, &report)
	if report.Equal {
		t.Fatalf("comparison report = %+v, want inequality", report)
	}
	if report.FirstDifference == "" {
		t.Fatal("comparison report omitted first difference")
	}
}

func writeFixtureBAM(t *testing.T, bamPath string) {
	t.Helper()
	header := &bamnative.Header{
		Version:   "1.6",
		SortOrder: "coordinate",
		References: []*bamnative.Reference{
			{ID: 0, Name: "chr1", Len: 1000},
		},
		OtherHeaderLines: []string{
			"@RG\tID:fixture\tSM:sample",
			"@PG\tID:fixture\tPN:bamroundtrip",
		},
	}
	writer, err := bamnative.NewWriter(bamPath, header)
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}
	fixtureRecords := []*bamnative.Record{
		{
			Name:      "fragment/1",
			Flags:     bamnative.FlagPaired | bamnative.FlagFirstInPair,
			RefID:     0,
			Pos:       10,
			MapQ:      60,
			Cigar:     []bamnative.CigarOp{{Op: bamnative.CigarSoftClip, Len: 1}, {Op: bamnative.CigarMatch, Len: 3}, {Op: bamnative.CigarInsertion, Len: 1}},
			MateRefID: 0,
			MatePos:   30,
			TLen:      25,
			Seq:       "TACGT",
			Qual:      []byte{30, 31, 32, 33, 34},
			Aux: []*bamnative.AuxField{
				{Tag: "NM", Type: bamnative.AuxTypeInt32, Value: int32(2)},
				{Tag: "XA", Type: bamnative.AuxTypeArray, ArrayType: bamnative.AuxTypeUInt8, Value: []byte{1, 2, 3}},
				{Tag: "ZF", Type: bamnative.AuxTypeFloat, Value: float32(1.5)},
			},
		},
		{
			Name:      "fragment/2",
			Flags:     bamnative.FlagPaired | bamnative.FlagSecondInPair | bamnative.FlagReverse,
			RefID:     0,
			Pos:       30,
			MapQ:      42,
			Cigar:     []bamnative.CigarOp{{Op: bamnative.CigarEqual, Len: 4}, {Op: bamnative.CigarDeletion, Len: 2}, {Op: bamnative.CigarMismatch, Len: 1}},
			MateRefID: 0,
			MatePos:   10,
			TLen:      -25,
			Seq:       "ACGTT",
			Aux: []*bamnative.AuxField{
				{Tag: "NM", Type: bamnative.AuxTypeInt32, Value: int32(3)},
				{Tag: "SA", Type: bamnative.AuxTypeString, Value: "chr1,50,+,5M,60,0;"},
			},
		},
	}
	for recordIndex, record := range fixtureRecords {
		if err := writer.Write(record); err != nil {
			_ = writer.Close()
			t.Fatalf("Write record %d: %v", recordIndex, err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

func writeMutatedFixtureBAM(t *testing.T, bamPath string) {
	t.Helper()
	writeFixtureBAM(t, bamPath)

	inputFile, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("open fixture: %v", err)
	}
	reader, err := bamnative.NewReader(inputFile)
	if err != nil {
		_ = inputFile.Close()
		t.Fatalf("NewReader: %v", err)
	}
	firstRecord, err := reader.Read()
	if err != nil {
		_ = inputFile.Close()
		t.Fatalf("read first record: %v", err)
	}
	secondRecord, err := reader.Read()
	if err != nil {
		_ = inputFile.Close()
		t.Fatalf("read second record: %v", err)
	}
	if err := inputFile.Close(); err != nil {
		t.Fatalf("close fixture: %v", err)
	}

	firstRecord.MapQ = 1
	writer, err := bamnative.NewWriter(bamPath, reader.Header())
	if err != nil {
		t.Fatalf("NewWriter mutation: %v", err)
	}
	if err := writer.Write(firstRecord); err != nil {
		_ = writer.Close()
		t.Fatalf("Write mutated record: %v", err)
	}
	if err := writer.Write(secondRecord); err != nil {
		_ = writer.Close()
		t.Fatalf("Write unmodified record: %v", err)
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close mutated fixture: %v", err)
	}
}

func readReport(t *testing.T, reportPath string, target interface{}) {
	t.Helper()
	reportBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("ReadFile report: %v", err)
	}
	if err := json.Unmarshal(reportBytes, target); err != nil {
		t.Fatalf("Unmarshal report: %v", err)
	}
}
