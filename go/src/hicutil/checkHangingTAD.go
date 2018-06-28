package hicutil

import (
//	"fmt"
//	"os"
)

func ContainsHangingTAD(tadlist [][]int, start int, end int) bool {

	for _,tad := range tadlist {
		if start > tad[1] {
			continue
		}
		if start >= tad[0] && start <= tad[1] { // interval ends in current TAD
			if start >= 3*(tad[1] - tad[0])/4 + tad[0] {
				return true
			}
		} else if end >= tad[0] && end <= tad[1] { // interval ends in current TAD
			if end <= 3*(tad[1] - tad[0])/4 + tad[0] {
				return true
			}
		} else if end <= tad[0] { // past end of interval
			break
		}

	}
	return false
}
