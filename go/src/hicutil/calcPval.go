package hicutil

import (
	"math"
	"math/big"
	"math/rand"
//	"sync"
	"runtime"
//	"fmt"
//	"os"
)

func CalcPval(intvl1 [][]int, intvl2 [][]int, n int, vival float64) float64 {
// more efficient to calculate VI as entropy(intvl1) + entropy(intvl2) - 2*mutinfo, because only need to recalculate mutual info on each iteration
	nshuffles := 10000
	if len(intvl1) < 8 && len(intvl2) < 8 {
		var f big.Int
		f.MulRange(1, int64(math.Max(float64(len(intvl1)), float64(len(intvl2)))))
		nshuffles = (int)(f.Int64())
	}
	//shuffvi1 := make([]float64, nshuffles)
	//shuffvi2 := make([]float64, nshuffles)

	h1 := CalcEntropy(intvl1,n)
	h2 := CalcEntropy(intvl2,n)

	clus1sizes := make([]int, len(intvl1))
        for c,clus := range intvl1 {
		clus1sizes[c] = clus[1]-clus[0]+1 }
	clus2sizes := make([]int, len(intvl2))
	for c,clus := range intvl2 {
		clus2sizes[c] = clus[1]-clus[0]+1 }

	count1 := 0
        count2 := 0
	for i := 0; i < nshuffles; i++ {
		runtime.Gosched()
		// randomly shuffle domain lengths in each list
		newlist1,permsizes1 := shuffledoms(clus1sizes)
		newlist2,permsizes2 := shuffledoms(clus2sizes)

		//calc VI for newly shuffled domains
		overlaps1 := CalcOverlaps(intvl1, newlist2)
		overlaps2 := CalcOverlaps(newlist1, intvl2)
		mutinfo1 := CalcMutInfo(overlaps1, clus1sizes, permsizes2, n)
		mutinfo2 := CalcMutInfo(overlaps2, permsizes1, clus2sizes, n)

		shuffvi1 := (h1+h2-2*mutinfo1)/math.Log(float64(n)) // divide by log(n) to normalize
		shuffvi2 := (h1+h2-2*mutinfo2)/math.Log(float64(n)) // divide by log(n) to normalize

		if shuffvi1 - vival < 1e-10 {
			count1++
		}
		if shuffvi2 - vival < 1e-10 {
                        count2++
                }

	}
	// find how many times shuffvi values are less than given vival
	/*count1 := 0
	count2 := 0
	for i:=0; i< nshuffles; i++ {
		if shuffvi1[i] - vival < 1e-10 { count1++ }
		if shuffvi2[i] - vival < 1e-10 { count2++ }
	}*/
	pval := (float64(count1+1)/float64(nshuffles+1) + float64(count2+1)/float64(nshuffles+1))/2.0
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
