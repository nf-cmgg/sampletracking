package main

import (
	"testing"
)

func TestFindFilesByExtension(t *testing.T) {
	dir := "test-data"
	extension := []string{".txt", ".csv"}

	files, err := findFilesByExtension(dir, extension)
	if err != nil {
		t.Errorf("Error finding files: %v", err)
	}

	expectedFiles := []string{
		"test-data/test1.txt",
		"test-data/test2.csv",
		// Add more expected file paths here
	}

	if len(files) != len(expectedFiles) {
		t.Errorf("Expected %d files, but got %d", len(expectedFiles), len(files))
	}

	for i, file := range files {
		if file != expectedFiles[i] {
			t.Errorf("Expected file path %s, but got %s", expectedFiles[i], file)
		}
	}
}
