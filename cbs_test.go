package cbsgo_test

import (
	"reflect"
	"sort"
	"testing"

	"github.com/matthds/cbsgo"
)

func TestCBS(t *testing.T) {
	// The input data from the original Rust test.
	// Since our function accepts float64, we convert the data.
	steps := []float64{1, 1, 1, 3, 3, 2, 1, 2, 3, 300, 310, 321, 310, 299}
	shuffles := 1000
	p := 0.05
	seed := int64(42) // Use a fixed seed for reproducibility.

	// Expected output. The segmentation may vary slightly due to randomness.
	// The main split is between the low and high part.
	expected := [][2]int{{0, 8}, {8, 14}}

	res, err := cbsgo.CBS(steps, shuffles, p, seed)
	if err != nil {
		t.Fatalf("CBS function returned an unexpected error: %v", err)
	}

	// Sort the resulting segments by their start point for consistent comparison.
	sort.Slice(res, func(i, j int) bool {
		return res[i][0] < res[j][0]
	})

	if !reflect.DeepEqual(res, expected) {
		t.Errorf("Unexpected result.\nExpected: %v\nGot: %v", expected, res)
	}
}
