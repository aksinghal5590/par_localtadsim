package hicutil

import (
	"math"
//	"math/big"
	"math/rand"
	"sync"
	"runtime"
//	"fmt"
//	"os"
)

func CalcPval(intvl1 [][]int, intvl2 [][]int, n int, vival float64, threadCount int) float64 {

	nshuffles := 10000/threadCount
	shuffvi1 := make([][]float64, threadCount)
	shuffvi2 := make([][]float64, threadCount)

	// more efficient to calculate VI as entropy(intvl1) + entropy(intvl2) - 2*mutinfo, because only need to recalculate mutual info on each iteration
	h1 := CalcEntropy(intvl1, n)
	h2 := CalcEntropy(intvl2, n)

	clus1sizes := make([]int, len(intvl1))
		for c,clus := range intvl1 {
		clus1sizes[c] = clus[1]-clus[0]+1
	}
	clus2sizes := make([]int, len(intvl2))
	for c,clus := range intvl2 {
		clus2sizes[c] = clus[1]-clus[0]+1
	}

	var wg sync.WaitGroup
	for w := 0; w < threadCount; w++ {
		wg.Add(1)
		shuffvi1[w] = make([]float64, nshuffles)
		shuffvi2[w] = make([]float64, nshuffles)
		go func(shuff1 []float64, shuff2 []float64) {
			for i := 0; i < nshuffles; i++ {
				runtime.Gosched()
				shuff1[i], shuff2[i] = worker(intvl1, intvl2, n, clus1sizes, clus2sizes, h1, h2)
			}
			defer wg.Done()
		}(shuffvi1[w], shuffvi2[w])
	}
	wg.Wait()

	// find how many times shuffvi values are less than given vival
	count1 := 0
	count2 := 0

	for w := 0; w < threadCount; w++ {
		for i := 0; i < nshuffles; i++ {
			if shuffvi1[w][i] - vival < 1e-10 {
				count1++
			}
			if shuffvi2[w][i] - vival < 1e-10 {
				count2++
			}
		}
	}
	pval := (float64(count1 + 1)/float64(nshuffles * threadCount + 1) + float64(count2 + 1)/float64(nshuffles * threadCount + 1))/2.0
	return pval
}

func worker(intvl1 [][]int, intvl2 [][]int, n int, clus1sizes []int, clus2sizes []int, h1 float64, h2 float64) (float64, float64) {

	// randomly shuffle domain lengths in each list
	newlist1, permsizes1 := shuffledoms(clus1sizes)
	newlist2, permsizes2 := shuffledoms(clus2sizes)

	//calc VI for newly shuffled domains
	overlaps1 := CalcOverlaps(intvl1, newlist2)
	overlaps2 := CalcOverlaps(newlist1, intvl2)
	mutinfo1 := CalcMutInfo(overlaps1, clus1sizes, permsizes2, n)
	mutinfo2 := CalcMutInfo(overlaps2, permsizes1, clus2sizes, n)

	shuffvi1 := (h1+h2-2*mutinfo1)/math.Log(float64(n)) // divide by log(n) to normalize
	shuffvi2 := (h1+h2-2*mutinfo2)/math.Log(float64(n)) // divide by log(n) to normalize

	return shuffvi1, shuffvi2
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
