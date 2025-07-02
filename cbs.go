package cbs

import (
	"math"
	"math/rand"
	"time"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/stat"
)

// CBS performs the Circular Binary Segmentation algorithm.
// It segments the input data `x` into pieces with a significantly different mean.
//
// Parameters:
//   - x: A slice of float64 data.
//   - shuffles: The number of permutations to perform to determine significance (1000 is recommended).
//   - significanceLevel: The p-value significance level (0.05 is recommended).
//
// Returns:
//   - A slice of [2]int arrays, where each array represents a [start, end] interval of a segment.
//   - An error if something goes wrong during the calculation.
func CBS(x []float64, shuffles int, significanceLevel float64, seed int64) ([][2]int, error) {
	// Use a seeded random source for reproducible shuffles.
	// For true randomness, use a different seed, e.g., time.Now().UnixNano().
	var rng *rand.Rand
	if seed == 0 {
		rng = rand.New(rand.NewSource(seed))
	} else {
		rng = rand.New(rand.NewSource(time.Now().UnixNano()))
	}

	var segments [][2]int
	err := rsegment(x, 0, len(x), &segments, shuffles, significanceLevel, rng)
	if err != nil {
		return nil, err
	}
	return segments, nil
}

// rsegment is the recursive function that performs the segmentation.
func rsegment(x []float64, start, end int, l *[][2]int, shuffles int, p float64, rng *rand.Rand) error {
	if start >= end {
		return nil
	}

	isChange, _, s, e, err := cbsInner(x[start:end], shuffles, p, rng)
	if err != nil {
		return err
	}

	// Add segment if there is no significant changepoint or if the segment is too small.
	if !isChange || (e-s < 5) || (e-s == end-start) {
		*l = append(*l, [2]int{start, end})
		return nil
	}

	// Recursively call for the sub-segments.
	// Segment before the changepoint
	if s > 0 {
		if err := rsegment(x, start, start+s, l, shuffles, p, rng); err != nil {
			return err
		}
	}
	// Segment of the changepoint itself
	if e-s > 0 {
		*l = append(*l, [2]int{start + s, start + e})
	}
	// Segment after the changepoint
	if start+e < end {
		if err := rsegment(x, start+e, end, l, shuffles, p, rng); err != nil {
			return err
		}
	}

	return nil
}

// cbsInner determines if there is a significant changepoint in the slice `x`.
func cbsInner(x []float64, shuffles int, p float64, rng *rand.Rand) (bool, float64, int, int, error) {
	maxT, maxStart, maxEnd, err := cbsStat(x)
	if err != nil {
		return false, 0, 0, 0, err
	}

	if maxEnd-maxStart == len(x) {
		return false, maxT, maxStart, maxEnd, nil
	}

	// Adjust start/end according to the heuristic in the original code.
	if maxStart < 5 {
		maxStart = 0
	}
	if len(x)-maxEnd < 5 {
		maxEnd = len(x)
	}

	// Permutation test
	threshCount := 0
	alpha := float64(shuffles) * p
	xt := make([]float64, len(x))
	copy(xt, x)

	for i := 0; i < shuffles; i++ {
		rng.Shuffle(len(xt), func(i, j int) { xt[i], xt[j] = xt[j], xt[i] })
		threshold, _, _, err := cbsStat(xt)
		if err != nil {
			return false, 0, 0, 0, err
		}
		if threshold >= maxT {
			threshCount++
		}
		if float64(threshCount) > alpha {
			return false, maxT, maxStart, maxEnd, nil
		}
	}

	return true, maxT, maxStart, maxEnd, nil
}

// cbsStat calculates the CBS test statistic.
// It uses gonum for efficient calculations.
func cbsStat(x []float64) (float64, int, int, error) {
	if len(x) == 0 {
		return 0.0, 0, 0, nil
	}

	length := float64(len(x))
	mean := stat.Mean(x, nil)

	// Create a mean-centered slice
	x0 := make([]float64, len(x))
	for i, val := range x {
		x0[i] = val - mean
	}

	// Calculate the cumulative sum with gonum
	y := make([]float64, len(x0))
	floats.CumSum(y, x0)

	// Find the indices of the max and min values in the cumulative sum
	e0 := floats.MaxIdx(y)
	e1 := floats.MinIdx(y)

	i0 := e0
	i1 := e1
	if e1 < e0 {
		i0 = e1
		i1 = e0
	}

	s0 := y[i0]
	s1 := y[i1]

	// Calculate the statistic
	denominator := (float64(i1-i0) + 1.0) * (length - float64(i1-i0))
	if denominator == 0 {
		return 0.0, i0, i1 + 1, nil // Avoid division by zero
	}
	stat := math.Pow(s1-s0, 2) * length / denominator

	return stat, i0, i1 + 1, nil
}
