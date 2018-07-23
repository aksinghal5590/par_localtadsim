package hicutil

import (
	"math"
	"math/rand"
)

func CalcPval(intvl1 [][]int, intvl2 [][]int, n int, vival float64, nshuffles int, r *rand.Rand) float64 {

	// more efficient to calculate VI as entropy(intvl1) + entropy(intvl2) - 2*mutinfo,
	// because only need to recalculate mutual info on each iteration
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

	overlaps1 := make([][]int, len(intvl1))
	for i,_ := range intvl1 {
		overlaps1[i] = make([]int, len(intvl2))
	}
	overlaps2 := make([][]int, len(intvl1))
	for i,_ := range intvl1 {
                overlaps2[i] = make([]int, len(intvl2))
        }

	count1 := 0
	count2 := 0

	shuffleCount := len(intvl1) * len(intvl2)
	if nshuffles > shuffleCount {
		nshuffles = shuffleCount
	}
	for i := 0; i < nshuffles; i++ {

		// randomly shuffle domain lengths in each list
		newlist1,permsizes1 := shuffledoms(clus1sizes, r)
		newlist2,permsizes2 := shuffledoms(clus2sizes, r)

		//calc VI for newly shuffled domains
		CalcOverlapsPtr(intvl1, newlist2, &overlaps1)
		CalcOverlapsPtr(newlist1, intvl2, &overlaps2)
		mutinfo1 := CalcMutInfo(overlaps1, clus1sizes, permsizes2, n)
		mutinfo2 := CalcMutInfo(overlaps2, permsizes1, clus2sizes, n)

		shuffvi1 := (h1+h2-2*mutinfo1)/math.Log(float64(n)) // divide by log(n) to normalize
		shuffvi2 := (h1+h2-2*mutinfo2)/math.Log(float64(n)) // divide by log(n) to normalize

		// find how many times shuffvi values are less than given vival
		if shuffvi1 - vival < 1e-10 {
			count1++
		}
		if shuffvi2 - vival < 1e-10 {
			count2++
		}
	}
	pval := (float64(count1+1)/float64(nshuffles+1) + float64(count2+1)/float64(nshuffles+1))/2.0
	return pval
}

func shuffledoms(clussizes []int, r *rand.Rand) ([][]int, []int) {

	permsizes := make([]int, len(clussizes))
	perm := r.Perm(len(clussizes))
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
