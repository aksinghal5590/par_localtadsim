package main

import (
	"hicutil"
	"flag"
	"strings"
	"fmt"
	"os"
	"bufio"
	"strconv"
	"sort"
	"math"
	"math/rand"
	"runtime"
	"runtime/pprof"
	"sync"
	"time"
	)

type bdyvi struct {
	start int
	end int
	vi float64
	pval float64
}

func sortAndPrint(bdyvis []bdyvi, isPrint bool) {
	sort.Slice(bdyvis, func(i,j int) bool {
		if bdyvis[i].start == bdyvis[j].start {return bdyvis[i].end < bdyvis[j].end} else {return bdyvis[i].start < bdyvis[j].start}
	})
	if isPrint {
		for i := 0; i < len(bdyvis); i++ {
			fmt.Printf("%d\t%d\t%6f\t%6f\n", bdyvis[i].start, bdyvis[i].end, bdyvis[i].vi, bdyvis[i].pval)
		}
	}
}


func main() {

	runtime.GOMAXPROCS(1)
	numCPU := runtime.NumCPU()
	// inputs should be 2 TAD files, a resolution parameter, and the output file name
	tadin := flag.String("tad", "", "comma-separated list of two TAD filenames or file patterns if optimizing gamma")
	gammain := flag.String("gamma", "", "if optimizing gamma, use 'opt,n' where n is the median TAD size to optimize for")
	res := flag.Int("res",1,"resolution of Hi-C data")
	outfile := flag.String("o","","output filename")
	pcount := flag.Int("p", numCPU, "no of processes")
	cpuprof := flag.String("cpu", "", "write cpu profile to file")

	flag.Parse()

	if *cpuprof != "" {
		f,_ := os.Create(*cpuprof)
			pprof.StartCPUProfile(f)
			defer pprof.StopCPUProfile()
	}

	tadfilelist := strings.Split(*tadin, ",")

	var gammaopt bool
	medtadlen := 0.0
	var err error
	if len(*gammain) > 0 {
		gammaopt = true
		gammadata := strings.Split(*gammain,",")
		medtadlen,err = strconv.ParseFloat(gammadata[1],64)
		if err != nil {
			fmt.Println("Error: couldn't convert median TAD length value to float, make sure to input i.e. '-gamma=opt,100' ")
			os.Exit(1)
		}
	} else {
		gammaopt = false
	}

	// read TAD files and process TAD lists to fill in non-TAD regions
	then := time.Now()
	tadlists := processTADLists(tadfilelist, res, gammaopt, medtadlen)
	duration := time.Since(then)
	fmt.Printf("processTADLists duration: %s\n", duration)
	// fmt.Println("done processing TAD lists, choosing optimal gamma")

	// calculate VI values at boundaries (using DP)
	then = time.Now()
	bdyvis := calcVIatBdys(tadlists)
	// bdyvis := calcVIatBdysNaive(tadlists)
	sortAndPrint(bdyvis, false)
	duration = time.Since(then)
	fmt.Printf("calcVIatBdys duration: %s\n", duration)

	// calculate all p-values, select significant points
	convergencecondition := 1e-5
	numCPU = *pcount
	runtime.GOMAXPROCS(numCPU)
	then = time.Now()
	sigpts := calcAllPvals(tadlists, bdyvis, numCPU, convergencecondition)
	duration = time.Since(then)
	fmt.Printf("calcAllPvals duration: %s\n", duration)
	// fmt.Println("done calculating all p-values")

	// identify dominating points from significant ones
	runtime.GOMAXPROCS(1)
	then = time.Now()
	dompts := findDomPts(sigpts)
	duration = time.Since(then)
	fmt.Printf("findDomPts duration: %s\n", duration)
	// fmt.Println("done finding dominating points")

	// save results to a file
	writeOutputToFile(dompts,outfile)
}


func processTADLists(tadfilelist []string, res *int, gammaopt bool, medtadlen float64) ([][][]int) {

	tadlists := make([][][]int, 2)
	chrlength := 0
	var gamma float64
	for i:=0; i < 2; i++ {
		if gammaopt == true {
			tadlists[i],gamma = hicutil.ChooseGamma(medtadlen, tadfilelist[i], *res)
			_ = gamma
		} else {
			tadlists[i] = hicutil.ReadTADFile(tadfilelist[i], *res)
		}
		n := tadlists[i][len(tadlists[i])-1][1]
		if chrlength < n+1 {
			chrlength = n+1
		}
	}
	for i:=0; i < 2; i++ {
		tadlists[i] = hicutil.FillinTADList(tadlists[i], chrlength)
	}

	return tadlists

}

type condentropy struct {
	vi float64
	condh1 float64
	condh2 float64
}

func calcVIatBdysNaive(tadlists [][][]int) ([]bdyvi) {
	var bdyvilist []bdyvi
	var newbdyvi bdyvi
	for _,tadlist := range tadlists {
		for i,tadstart := range tadlist {
			for _,tadend := range tadlist[i:] {

				if hicutil.ContainsHangingTAD(tadlists[0], tadstart[0], tadend[1]) || hicutil.ContainsHangingTAD(tadlists[1], tadstart[0], tadend[1]) {
					continue
				}

				newbdyvi.start = tadstart[0]
				newbdyvi.end = tadend[1]
				n := tadend[1] - tadstart[0] + 1
				intvl1 := hicutil.ProcessIntervals(tadlists[0],tadstart[0],tadend[1])
				intvl2 := hicutil.ProcessIntervals(tadlists[1],tadstart[0],tadend[1])
				overlaps := hicutil.CalcOverlaps(intvl1,intvl2)
				clus1sizes := make([]int, len(intvl1))
				for c,clus := range intvl1 {
					clus1sizes[c] = clus[1]-clus[0]+1 }
				clus2sizes := make([]int, len(intvl2))
				for c,clus := range intvl2 {
					clus2sizes[c] = clus[1]-clus[0]+1 }
				condh1 := hicutil.CalcCondEntropy(transpose(overlaps), clus1sizes, n)
				condh2 := hicutil.CalcCondEntropy(overlaps, clus2sizes, n)
				newbdyvi.vi = (condh1 + condh2)/math.Log(float64(n))
				bdyvilist = append(bdyvilist, newbdyvi)
			}
		}
	}
	return bdyvilist
}

func scaleLastOverlap(overlaps [][]int, maxn int, lasttad []int) [][]int {

	scaleby := float64(maxn - lasttad[1] + 1)/float64(maxn - lasttad[0] + 1)
	overlaps[len(overlaps)-1][len(overlaps[0])-1] = int(float64(overlaps[len(overlaps)-1][len(overlaps[0])-1]) * scaleby)
	return overlaps
}

func calcVIatBdys(tadlists [][][]int) ([]bdyvi) {

	var bdyvilist []bdyvi
	dpmap := make(map[int]map[int]condentropy) // keys should be [start][end]
	var currhvals condentropy
	// first initialize VI values of all single-TAD intervals
	for _,tadlist := range tadlists {
		for _,tads := range tadlist {

			var newbdyvi bdyvi
			newbdyvi.start = tads[0]
			newbdyvi.end = tads[1]
			//if tads[0] == tads[1] { continue }
			n := tads[1] - tads[0] + 1
			//newend := extendn(tadlists[int(math.Abs(float64(i-1)))], tads[1])
			//n := newend - tads[0] + 1
			intvl1 := hicutil.ProcessIntervals(tadlists[0], tads[0], tads[1])
			intvl2 := hicutil.ProcessIntervals(tadlists[1], tads[0], tads[1])
			// need cluster sizes and overlaps for conditional entropy calculations
			overlaps := hicutil.CalcOverlaps(intvl1, intvl2)
			clus1sizes := make([]int, len(intvl1))
			for c,clus := range intvl1 {
				clus1sizes[c] = clus[1]-clus[0]+1 }
			clus2sizes := make([]int, len(intvl2))
			for c,clus := range intvl2 {
				clus2sizes[c] = clus[1]-clus[0]+1 }
			currhvals.condh1 = hicutil.CalcCondEntropy(transpose(overlaps), clus1sizes, n)
			currhvals.condh2 = hicutil.CalcCondEntropy(overlaps, clus2sizes, n)
			currhvals.vi = currhvals.condh1 + currhvals.condh2
			if _, ok := dpmap[tads[0]]; !ok {
				dpmap[tads[0]] = make(map[int]condentropy)
			}
			dpmap[tads[0]][tads[1]] = currhvals
			if math.IsNaN(currhvals.vi) {
				fmt.Println("VI value of NaN in initialization")
				fmt.Println(currhvals)
				fmt.Println(tads)
				fmt.Println(intvl1)
				fmt.Println(tadlists[1])
				os.Exit(1)
			}
			newbdyvi.vi = currhvals.vi/math.Log(float64(tads[1]-tads[0]+1)) // divide by log(n) to normalize
			if tads[0] == tads[1] { continue }
			bdyvilist = append(bdyvilist, newbdyvi)
		}
	}
	for i,tadlist := range tadlists {
		var noti int
		if i == 0 { noti = 1 } else { noti = 0 }
		for numtads := 1; numtads < len(tadlist); numtads++ {
			for tadidx, starttad := range tadlist[:len(tadlist)-numtads] {
				endtad := tadlist[tadidx+numtads]
				n := endtad[1] - starttad[0] + 1
				prevn := endtad[0] - starttad[0]
				endtadlen := endtad[1]-endtad[0]+1
				scale1 := float64(prevn)/float64(n)
				scale2 := float64(endtadlen)/float64(n)
				var newhvals condentropy
				var newbdyvis bdyvi
				newbdyvis.start = starttad[0]
				newbdyvis.end = endtad[1]
				if contains(tadlists[noti], endtad[0]-1) {
					// there are no overlapping TADs between the previous group and the new TAD so can just rescale and add the VIs/conditional entropies
					newhvals.vi =  dpmap[starttad[0]][endtad[0]-1].vi*scale1 + dpmap[endtad[0]][endtad[1]].vi*scale2
					newhvals.condh1 = dpmap[starttad[0]][endtad[0]-1].condh1*scale1+dpmap[endtad[0]][endtad[1]].condh1*scale2
					newhvals.condh2 = dpmap[starttad[0]][endtad[0]-1].condh2*scale1+dpmap[endtad[0]][endtad[1]].condh2*scale2
				} else { // the other TAD set contains a TAD that crosses the boundary between endtad-1 and endtad, so we have to correct one of the conditional entropy terms
					if i == 0 {
						newhvals.condh1 = scale1*dpmap[starttad[0]][endtad[0]-1].condh1 + scale2*dpmap[endtad[0]][endtad[1]].condh1
						newhvals.condh2 = CorrectTADoverlap(tadlists, i, starttad[0], endtad, dpmap)
					} else {
						newhvals.condh2 = scale1*dpmap[starttad[0]][endtad[0]-1].condh2 + scale2*dpmap[endtad[0]][endtad[1]].condh2
						newhvals.condh1 = CorrectTADoverlap(tadlists, i, starttad[0], endtad, dpmap)
					}
					newhvals.vi = newhvals.condh1 + newhvals.condh2
				}
				if math.IsNaN(newhvals.vi) {
					fmt.Println("VI value of NaN in DP")
					fmt.Println(newhvals)
					fmt.Println(starttad, endtad)
					os.Exit(1)
				}
				dpmap[starttad[0]][endtad[1]] = newhvals
				newbdyvis.vi = newhvals.vi/math.Log(float64(endtad[1]-starttad[0]+1))
				if hicutil.ContainsHangingTAD(tadlists[0], starttad[0], endtad[1]) || hicutil.ContainsHangingTAD(tadlists[1], starttad[0], endtad[1]) {
					continue
				}
				bdyvilist = append(bdyvilist, newbdyvis)
			}
		}
	}
	return bdyvilist
}

func contains(a [][]int, x int) bool {
	for _, row := range a {
		if row[1] == x {
			return true
		}
	}
	return false
}

func transpose(a [][]int) [][]int {
	n := len(a)
	m := len(a[0])
	b := make([][]int, m)
	for j := 0; j < m; j++ {
		b[j] = make([]int, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			b[j][i] = a[i][j]
		}
	}
	return b
}

var wg sync.WaitGroup
func worker(tadlists [][][]int, job []bdyvi, result *[]bdyvi, convcond float64) {
	defer wg.Done()
	concurrentRand := rand.New(rand.NewSource(time.Now().UnixNano()))
	for i, querypt := range job {
		(*result)[i] = appendPval(tadlists, querypt, convcond, concurrentRand)
	}
}

func calcAllPvals(tadlists [][][]int, bdyvis []bdyvi, numCPU int, convcond float64) []bdyvi {

	var sigpts []bdyvi
	if len(bdyvis) < numCPU {
		numCPU = len(bdyvis)
	}
	bdyvis_pval := make([]bdyvi, len(bdyvis))
	allpvals := make([]float64, len(bdyvis))
	jobs := make([][]bdyvi, numCPU)
	results := make([][]bdyvi, numCPU)

	k := 0
	for i := 0; i < len(bdyvis); i++ {
		jobs[k] = append(jobs[k], bdyvis[i])
		k = (k + 1) % numCPU
	}

	for w, job := range jobs {
		wg.Add(1)
		results[w] = make([]bdyvi, len(job))
		go worker(tadlists, job, &results[w], convcond)
	}
	wg.Wait()

	i := 0
	for w := 0; w < numCPU; w++ {
		for _, res := range results[w] {
			bdyvis_pval[i] = res
			allpvals[i] = bdyvis_pval[i].pval
			i++
		}
	}

	bhidx := hicutil.MultHypTestBH(allpvals)
	if bhidx > -1 {
		sort.Slice(bdyvis_pval, func(i,j int) bool {return bdyvis_pval[i].pval < bdyvis_pval[j].pval})
		sigpts = bdyvis_pval[:bhidx+1]
	}
	return sigpts
}

/*func calcAllPvals(tadlists [][][]int, bdyvis []bdyvi, nshuffles int) []bdyvi {

	var sigpts []bdyvi
	bdyvis_pval := make([]bdyvi, len(bdyvis))
	allpvals := make([]float64,len(bdyvis_pval))
	for i,querypt := range bdyvis {
		bdyvis_pval[i] = appendPval(tadlists, querypt, nshuffles)
		allpvals[i] = bdyvis_pval[i].pval
	}
	bhidx := hicutil.MultHypTestBH(allpvals)
	sort.Slice(bdyvis_pval, func(i,j int) bool {return bdyvis_pval[i].pval < bdyvis_pval[j].pval})
	sigpts = bdyvis_pval[:bhidx+1]
	return sigpts
}*/

func appendPval( tadlists [][][]int, querypt bdyvi, convcond float64, r *rand.Rand) (bdyvi) {

	intvl1 := hicutil.ProcessIntervals(tadlists[0], querypt.start, querypt.end)
	intvl2 := hicutil.ProcessIntervals(tadlists[1], querypt.start, querypt.end)
	n := querypt.end - querypt.start + 1
	p := hicutil.CalcPval(intvl1, intvl2, n, querypt.vi, convcond, r)
	if p < 0.05 && (len(intvl1) == 1 || len(intvl2) == 1) {
		fmt.Println(intvl1)
		fmt.Println(intvl2)
		fmt.Println(p)
		fmt.Println(n,querypt.vi)
		fmt.Println(querypt.start, querypt.end)
		os.Exit(1)
	}
	querypt.pval = p
	return querypt
}

func findDomPts(sigpts []bdyvi) []bdyvi {

	var dompts []bdyvi
	for i,querypt := range sigpts {
		isdom := true
		for j,comppt := range sigpts {
			if i == j { continue }
			if comppt.start == querypt.start && comppt.end > querypt.end && comppt.vi < querypt.vi {
				isdom = false
				break
			}
			if comppt.start < querypt.start && comppt.end == querypt.end && comppt.vi < querypt.vi {
				isdom = false
				break
			}
			if comppt.start <= querypt.start { continue }
			if comppt.end >= querypt.end { continue }
			if comppt.vi < querypt.vi {
				isdom = false
				break
			}
		}
		if isdom == true {
			dompts = append(dompts, querypt)
		}
	}
	//if there are multiple dominating, significant points that start or end at the same place, remove the smaller interval
	var toremove []int
	for i,dompt1 := range dompts {
		for j,dompt2 := range dompts {
			if i == j {continue}
			if dompt1.start == dompt2.start && dompt1.end < dompt2.end {
				toremove = append(toremove, i)
				break
			}
			if dompt1.start > dompt2.start && dompt1.end == dompt2.end {
				toremove = append(toremove, i)
				break
			}
			if dompt1.start == dompt2.start && dompt1.end == dompt2.end && i < j {
				toremove = append(toremove, i)
				break
			}
		}
	}
	sort.Ints(toremove)
	for i,j := range toremove {
		dompts = append(dompts[:j-i], dompts[j-i+1:]...)
	}
	sort.Slice(dompts, func(i,j int) bool {return dompts[i].start < dompts[j].start})
	return dompts
}

func writeOutputToFile(domsigpts []bdyvi, outfile *string) {

    //write values to file
    f,err := os.Create(*outfile)
    if err != nil {
            panic(err)
    }
    //defer f.Close()

    w := bufio.NewWriter(f)
    labelline := []string{"start", "end", "VI", "p-value"}
    //fmt.Println(strings.Join(labelline, "\t"))
	line1 := strings.Join(labelline, "\t")
    //fmt.Println(line1)
    fmt.Fprintf(w,line1+"\n")

    for _,vals := range domsigpts {
            strvals := make([]string, 4)
	strvals[0] = strconv.Itoa(vals.start)
	strvals[1] = strconv.Itoa(vals.end)
	strvals[2] = strconv.FormatFloat(vals.vi,'g',6,64)
	strvals[3] = strconv.FormatFloat(vals.pval,'g',6,64)
            newline := strings.Join(strvals, "\t")
            fmt.Fprintf(w,newline+"\n")
    }
    w.Flush()
    f.Close()
	//fmt.Println("Wrote output values to", *outfile)
}

