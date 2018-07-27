package hicutil

import (
	"math"
	"math/rand"
//	"sync"
//	"runtime"
//	"fmt"
//	"os"
)

func CalcPval(intvl1 [][]int, intvl2 [][]int, n int, vival float64, convcond float64) float64 {
// more efficient to calculate VI as entropy(intvl1) + entropy(intvl2) - 2*mutinfo, because only need to recalculate mutual info on each iteration
	h1 := CalcEntropy(intvl1,n)
	h2 := CalcEntropy(intvl2,n)

	clus1sizes := make([]int, len(intvl1))
        for c,clus := range intvl1 {
		clus1sizes[c] = clus[1]-clus[0]+1 }
	clus2sizes := make([]int, len(intvl2))
	for c,clus := range intvl2 {
		clus2sizes[c] = clus[1]-clus[0]+1 }

	var pval1 float64
	var pval2 float64
	count := 0
	shuffnum := 0
	prevpvaldiffs := []float64{1000.0, 1000.0, 1000.0, 1000.0}
	prevpval := 1000.0
	keepshuffling := true

	for keepshuffling {
		shuffnum++

		newlist, permsizes := shuffledoms(clus1sizes)
		overlaps := CalcOverlaps(newlist, intvl2)
		mutinfo := CalcMutInfo(overlaps, permsizes, clus2sizes, n)
		shuffvi := (h1+h2-2*mutinfo)/math.Log(float64(n)) // divide by log(n) to normalize

		// calc current pval, and pvaldiffs
		if shuffvi - vival < 1e-10 {
			count++
		}
		pval1 = float64(count+1)/float64(shuffnum+1)
		prevpvaldiffs = append(prevpvaldiffs[1:], []float64{math.Abs(prevpval-pval1)}...)

		// if last 5 p-values are within convergence condition, we're done shuffling
		keepshuffling = false
		for _,pvaldiff := range prevpvaldiffs {
			if pvaldiff > convcond {
				keepshuffling = true
			}
		}
		prevpval = pval1
	}

	// re-initialize variables
	shuffnum = 0
	count = 0
	prevpvaldiffs = []float64{1000.0, 1000.0, 1000.0, 1000.0}
	prevpval = 1000.0
	keepshuffling = true

	for keepshuffling {
		shuffnum++

		newlist, permsizes := shuffledoms(clus2sizes)
		overlaps := CalcOverlaps(intvl1, newlist)
		mutinfo := CalcMutInfo(overlaps, clus1sizes, permsizes, n)
		shuffvi := (h1+h2-2*mutinfo)/math.Log(float64(n)) // divide by log(n) to normalize

		// calc current pval, and pvaldiffs
		if shuffvi - vival < 1e-10 {
			count++
		}
		pval2 = float64(count+1)/float64(shuffnum+1)
		prevpvaldiffs = append(prevpvaldiffs[1:], []float64{math.Abs(prevpval-pval2)}...)

		// if last 5 p-values are within convergence condition, we're done shuffling
		keepshuffling = false
		for _,pvaldiff := range prevpvaldiffs {
			if pvaldiff > convcond {
				keepshuffling = true
			}
		}
		prevpval = pval2
	}

	pval := (pval1 + pval2)/2.0
	return pval
}

func shuffledoms(clussizes []int) ([][]int, []int) {

	permsizes := make([]int, len(clussizes))
	perm := rand.Perm(len(clussizes))
	for j,v := range perm {
		permsizes[v] = clussizes[j]
	}
	// turn shuffled lists of lengths back into TAD lists
	newlist := make([][]int, len(clussizes))
	newlist[0] = []int{0,permsizes[0]-1}
	for j := 1; j < len(permsizes); j++ {
		tadstart := newlist[j-1][1]+1
		newlist[j] = []int{tadstart, tadstart+permsizes[j]-1}
	}

	return newlist, permsizes
}
