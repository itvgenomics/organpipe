"""
@author: maze
"""

import copy
import glob
import logging
import math
import os
import sys
from typing import Dict, List, Literal

import Bio.Blast.NCBIXML
from Bio.Data import CodonTable

from mitos import feature, trna, update
from mitos.mergefeature import _copy_part
from mitos.mito import type_from_name
from mitos.sequence import sequences_fromfile


class UpdateZeroScore(Exception):
    """
    Exception to be raised when no _update_hit results in 0 score.
    """

    def __init__(self, hit):
        self._hit = hit

    def __str__(self):
        """
        print exception
        """
        return "_update_hit resulted in 0 score for %s" % str(self._hit)


class position_values:
    """
    data structure to store the values to be stored for a position. i.e.
    - its position (nucleotide position)
    - height
    - sum of evalues
    - sum of bitscores
    - query relative query position
    - yreuq relative query position relative to end
    """

    def __init__(self, pos):
        """
        initialise everything to 0 and the position to the given value
        @param pos: the position
        """
        self._position = pos
        self._height = 0
        self._evalue = 0.0
        self._bitscore = 0.0
        #         self._query = decimal.Decimal( 0.0 )
        self._query = 0.0
        self._yreuq = 0.0

    #         self._query = fractions.Fraction( 0, 1 )
    #         self._query = decimal.Decimal( 0 )

    def __str__(self):
        return "pv(%d %d %.1f %f %f)" % (
            self._position,
            self._height,
            round(self._evalue),
            self._query,
            self._yreuq,
        )

    def __repr__(self):
        return "pv(%d %d %.1f %f %f)" % (
            self._position,
            self._height,
            round(self._evalue),
            self._query,
            self._yreuq,
        )

    def add(self, evalue, bitscore, relquery, relyreuq):
        """
        add values of a given blast hit and increase height by one
        @param evalue: the evalue of the blast hit
        @param bitscore:  the evalue of the blast hit
        @param relquery: the relative query of the position in the blast hit
        @param relyreuq: the relative query of the position in the blast hit (wrt end)
        """

        self._height += 1
        self._evalue += evalue
        self._bitscore += bitscore
        self._query += relquery
        self._yreuq += relyreuq

    def get_bitscore(self, avg):
        """
        get the (average) bitscore
        @param avg: if true return avg, else return sum
        @return the desired value
        """
        if avg:
            try:
                return float(self._bitscore) / self._height
            except ZeroDivisionError:
                return 0
        else:
            return self._bitscore

    def get_height(self):
        """
        get the height
        @return the desired value
        """
        return self._height

    def get_evalue(self, avg):
        """
        get the (average) evalue
        @param avg: if true return avg, else return sum
        @return the desired value
        """
        if avg:
            try:
                return self._evalue / self._height
            except ZeroDivisionError:
                return 0
        else:
            return self._evalue

    def get_position(self):
        """
        just get the position
        @return: the position
        """
        return self._position

    def get_query(self):
        """
        get the average query position
        @return the desired value
        """
        return self._query

    def get_score(self, sel="e", avg=False):
        """
        get one of the three values
        @param sel: selector 'e': get evalue, 'b': get bitscore, 'h': get height
        @param avg: get the average value (useless for height)
        """
        if sel == "e":
            return self.get_evalue(avg)
        elif sel == "b":
            return self.get_bitscore(avg)
        elif sel == "h":
            return self.get_height()
        else:
            raise Exception("unknown selector" + repr(sel))

    def get_yreuq(self):
        """
        get the average query position (wrt end)
        @return the desired value
        """
        return self._yreuq

    def finalize(self):
        """
        finalize the position, i.e. average the query
        """
        self._query /= float(self._height)
        self._yreuq /= float(self._height)

    def set_position(self, p):
        """
        just set the position
        @param[in] p position
        @return: the position
        """
        self._position = p

    def set_query(self, q):
        """
        just set the query
        @param[in] p position
        @return: the position
        """
        self._query = q

    def sub(self, evalue, bitscore, relquery, relyreuq):
        """
        subtract values of a given blast hit and increase height by one
        @param evalue: the evalue of the blast hit
        @param bitscore:  the evalue of the blast hit
        @param relquery: the relative query of the position in the blast hit
        @param relquery: the relative query of the position in the blast hit (wrt end)
        """

        self._height -= 1
        self._evalue -= evalue
        self._bitscore -= bitscore
        self._query -= relquery
        self._yreuq -= relyreuq

    def set_yreuq(self, q):
        """
        just set the query
        @param[in] p position
        @return: the position
        """
        self._yreuq = q


def _apply_cutoff(idx, values, cutoff):
    """
    apply cutoff in order to make the hills steep.
    remove positions with less than (maxhits in the reading frame)*cutoff
    from the index and invalidate the corresponding entries in the
    three arrays (not deleted there because the mapping would become invalid)
    currently the cutoff is computed relative to the maximum number of hits in the rf

    TODO(later) do this per 'hit-component' otherwise pseudo genes with to few hits
    are removed

    @param idx an mapping from positions in the genome to an index for one of the genes
    @param values list of values per position (access via idx)
    @param cutoff cutoff cutoff value

    @return Nothing, idx and values is modified
    """
    for rf in list(idx.keys()):
        mx = 0
        for i in list(idx[rf].keys()):
            mx = max(mx, values[idx[rf][i]].get_score())
        endcut = mx * cutoff

        for i in list(idx[rf].keys()):
            p = idx[rf][i]
            if values[p].get_score() < endcut:
                values[p] = None
                del idx[rf][i]
        #            else:
        #                print p, i, values[p].get_score( )

        if len(idx[rf]) == 0:
            del idx[rf]


def blast(
    fastafiles,
    offsets,
    prot,
    cr,
    filepath,
    refdir,
    length,
    startcod,
    stopcod,
    fas2len,
    len_pval,
    idx,
    values,
    fragfac,
    fragovl,
    trnas,
    circular,
    cutoff=0,
    minevalue=0,
    acc=None,
    code=None,
    plot=False,
    prntih=False,
    finovl=100,
    clipfac=10,
    oldstst=False,
):
    """
    wrapper function to call
    1. singleblastx
    2. singleblastn
    3. blastx

    @param fastafiles list of fastafiles, 1 per offset
    @param offsets list of offsets for the fastafiles
    @param prot search for protein features
    @param cr search for nucleotide features
    @param filepath
    @param refdir
    @param length sequence length
    @param startcod start codon statistics
    @param stopcod stop codon statistics
    @param[out] idx
    @param[out] values
    @param featlength dictionary containing the lengths of the query sequences in the used DB
    @param len_pval a dictionary containing for each gene the length pvalues
    @param[in] oldstst use the old start stop prediction method of MITOS1
    """

    #        , ststrange = 6, maxdegenerate = 2,
    features = []
    if len(fastafiles) != len(offsets):
        raise Exception("blast.blast called with |sequences|!=|offsets|")

    # clear stst.dat
    open("%s/stst.dat" % (filepath), "w").close()

    for i in range(len(offsets) - 1, -1, -1):
        oidx = dict()
        oval = []

        if prot:
            update.singleblastx(fastafiles[i], code, filepath, refdir)
        if cr:
            update.singleblastn(fastafiles[i], filepath, refdir)

        # get an offset version of the forbidden regions
        otrnas = copy.deepcopy(trnas)
        for t in otrnas:
            t.start = (t.start - offsets[i]) % length
            t.stop = (t.stop - offsets[i]) % length

        # get features from blast results
        # - generate plots only for 1st offset
        tfl = blastx(
            "%sblast/" % (filepath),
            length,
            startcod,
            stopcod,
            fas2len,
            len_pval,
            oidx,
            oval,
            fragfac,
            fragovl,
            otrnas,
            circular=False,
            cutoff=cutoff,
            minevalue=minevalue,
            acc=acc,
            code=code,
            fastafile=fastafiles[i],
            plot=(i == 0),
            finovl=finovl,
            clipfac=clipfac,
            prot=prot,
            cr=cr,
            oldstst=oldstst,
        )
        #        ststrange = ststrange,maxdegenerate = maxdegenerate
        # correct the feature positions for the offset
        for f in tfl:
            f.start = (f.start + offsets[i]) % length
            f.stop = (f.stop + offsets[i]) % length

            f.rf = int(
                math.copysign(((abs(f.rf) - 1) + (offsets[i] % 3)) % 3 + 1, f.rf)
            )

            # the positions of the feature need not to be updated. this is done
            # below via values and the idx
            # for p in f._positions:
            #       p._position = ( p._position + offsets[i] ) % length
            features.append(f)
            logging.debug(f"offset {offsets[i]} {f}")

        # correct the idx and the values for the offset
        for n in oidx:  # names
            if n not in idx:
                idx[n] = dict()

            for of in oidx[n]:  # frames
                if abs(of) == 4:
                    f = of
                else:
                    f = int(
                        math.copysign(((abs(of) - 1) + (offsets[i] % 3)) % 3 + 1, of)
                    )

                if f not in idx[n]:
                    idx[n][f] = dict()

                for op in oidx[n][of]:  # positions
                    p = (op + offsets[i]) % length
                    v = oval[oidx[n][of][op]]
                    v.set_position((v.get_position() + offsets[i]) % length)

                    if p in idx[n][f]:
                        if v.get_score() > values[idx[n][f][p]].get_score():
                            values[idx[n][f][p]] = v
                    else:
                        idx[n][f][p] = len(values)
                        values.append(v)

    #                 logging.debug( "%s %d %s" % ( n, f, oidx[n][f].keys() ) )

    # output for plot histograms
    if plot:
        _create_protplotdata(filepath + "/blast/blast.dat", idx, values)

    return features


def blastx(
    filepath,
    length,
    startcod,
    stopcod,
    fas2len,
    len_pval,
    idx,
    values,
    fragfac,
    fragovl,
    trnas,
    circular,
    cutoff=0,
    minevalue=0,
    acc=None,
    code=None,
    fastafile=None,
    plot=False,
    prntih=False,
    finovl=100,
    clipfac=10,
    prot=True,
    cr=False,
    oldstst=False,
):
    """
    TODO description

    @param filepath path containing the blast results
        (prot/gen_name.blast and nuc/gen_name.blast files)
    @param length length of the sequence
    @param startcod start codon statistics
    @param stopcod stop codon statistics
    @param fas2len dictionary containing the lengths of the query sequences in the used DB
    @param len_pval a dictionary containing for each code and gene the
        mean (mu), standard deviation (sigma), maximum (max), number of data points (n) and pvalue list (pv)
    @param[out] idx an mapping from positions in the genome to an index in the following
        array. seperately for each gene and reading frame. idx[name][rf][i]
        gives the position in the arrays for position i for hits of gene with
        name in the specified reading frame
    @param[out] values the values of the blast search for each position
    @param cutoff the fraction [0:1] of quality value necessary to include
        a position in an initial hit
    @param minevalue
    @param acc accession/name of species (only used for determining plot path)
    @param code
    @param fastafile a the name of the file where the sequence can be found
           that generated the BLAST results
           (note this is the doubled sequence in case circular = True)
    @param plot
    @param prntih: print the names of the initial hits and their score and
           return afterwards
    @param clipfac of overlapping prediction of the same name differ by less
           than factor X than clipping is started
    @param prot predict proteins based on blastx results
    @param cr predict control region based on blastn
    @param[in] oldstst use the old start stop prediction method of MITOS1
    @return
    """
    # ststrange = 6, , maxdegenerate = 2
    #     @param ststrange # aa that should be searched up/downstream for a more precise
    #            start or stop position
    #     @param maxdegenerate maximum number of degenerated positions in stop codons

    #    print cutoff, minevalue, maxovl, clipfac, fragovl, fragfac, ststrange
    # print filepath, cutoff, minevalue, acc, code, fastafile, circular, plot,
    # debug, scoresel, pavg, havg, prntih

    # TODO remove acc (only used for determining the plot path -> store plots
    # in blast dir? )
    if acc is None:
        acc = fastafile.split("/")[-1]

    # if the code is given the algorithm searches for nearby start and stop codons
    # thus, the sequence has to be read in this case
    if code is not None:
        # read sequence
        sequence = sequences_fromfile(fastafile, circular=circular)[0]

    # list of the blast table files
    files = glob.glob("%s/*/%s.*.blast" % (filepath, os.path.basename(fastafile)))
    for f in files:
        # get name of gene (file name is xxx.GEN.blast)
        name = f.split(".")[-2]

        # determine if its in prot or nuc
        pth = f.split("/")[-2]

        if (pth == "nuc" and not cr) or (pth == "prot" and not prot):
            continue

        # read blast output
        idx[name] = _blasthandle(
            f, minevalue, values, fas2len, prot=(pth == "prot"), ftype="tab"
        )

    #         _apply_cutoff( idx[name], values, cutoff )
    #         break
    # END OF ITERATION OVER ALL blast files

    hits = _create_initial_predictions(idx, values, cutoff)

    logging.debug("initial hits")
    for h in hits:
        logging.debug("%s" % (str(h)))
    logging.debug("============")

    if prntih:
        for h in hits:
            sys.stdout.write("%s %f\n" % (h.name, h.score))

    # - 1st round of copy_part prediction (necessary because the start/stop
    # prediction is only for the first and last parts, i.e. this needs to be
    # known)
    cpyprt = _copy_part(hits, fragovl, fragfac)
    # - apply start and stop prediction to the 1st and last parts
    f = open("%s/stst.dat" % (filepath), "a")
    logging.debug("%s/stst.dat" % (filepath))
    for n in cpyprt:
        for cpy in cpyprt[n]:
            if cpyprt[n][cpy][0].type != "gene":
                continue
            if oldstst:
                _improve_start_stop_old(
                    cpyprt[n][cpy], code, sequence, len(sequence), idx[n], values
                )
            else:
                plotdat = _improve_start_stop(
                    cpyprt[n][cpy],
                    sequence,
                    code,
                    startcod,
                    stopcod,
                    len_pval,
                    idx[n],
                    values,
                    cpyprt,
                    trnas,
                    circular,
                )
                logging.debug("writing %d entries" % len(plotdat))
                f.write(plotdat)
    f.close()

    hits = _nonoverlapping_greedy(hits, idx, values, finovl, clipfac)
    #    hits = _nonoverlapping_greedy_sorted( hits, debug )

    # - run 2nd round of copy_part prediction, necessary because some parts will
    #   be removed by _nonoverlapping_greedy
    _copy_part(hits, fragovl, fragfac)

    #     # determine the best of each kind and store this info in the copy information
    #     best = {}
    #     for h in hits:
    #         if not h.name in best or h.score > best[h.name]:
    #             best[h.name] = h.score
    #
    #     for h in hits:
    #         if h.score == best[h.name]:
    #             h.copy = 0
    #         else:
    #             h.copy = 1

    return hits


def blastn(filepath, cutoff=0, minevalue=0, acc=None):
    if acc is None:
        acc = filepath.split("/")[-1]

    # Check if it is in %
    percent = False
    try:
        if cutoff[0] == "%":
            percent = True
            cutoff = float(cutoff[1:]) / 100
    except Exception:
        # if not percent the Cutoff ist a Int
        cutoff = int(cutoff)

    # list of the Posibel Gene
    posgenes = dict()
    strand = dict()
    # List of the Blasttable-Files
    protfiles = glob.glob("%s/*.blast" % filepath)
    # Overall blast Tables
    for _file in protfiles:
        try:
            # Get name of Gene
            name = _file.split(".")[-2]
            # Read file
            blasthandle = open(_file, "r")
            # Get Hits form the Table
            lines = blasthandle.readlines()
            # close File
            blasthandle.close()
            # Preper ReadingFreams
            posgenes[name] = {-1: dict(), 1: dict()}
            hits = {-1: 0, 1: 0}
            # For every Blast-Hit from the Gene
            for line in lines:
                # Splite the Table in the single Information
                cols = line.split()

                start = int(cols[6])
                stop = int(cols[7])
                # gettin evalue for the cut of evalue
                evalue = cols[10].split("-")[-1]
                if evalue == "0.0":
                    evalue = 100
                else:
                    evalue = float(evalue)

                # get strand
                if start < stop:
                    strand[name] = 1
                elif start > stop:
                    helpstart = start
                    start = stop
                    stop = helpstart
                    strand[name] = -1
                else:
                    logging.warning("[Blast] error: start = stop")
                    sys.exit()

                #####
                # check whether hit is identity, if so continue with next hit
                #####
                _id = False

                # identical accession number
                if cols[0] == cols[1].split(":")[0]:
                    _id = True
                    if not cols[2] == "100.00":
                        idstart = int(cols[1].split(":")[1].split("-")[0]) + 1
                        idstop = int(cols[1].split(":")[1].split("-")[1])
                        if idstop < idstart:
                            helpstop = idstart
                            idstart = idstop
                            idstop = helpstop
                        x = y = 0
                        for i in range(idstart, idstop):
                            y += 1
                            if i > start and i < stop:
                                x += 1  # count same positions
                        if float(x) / y < 0.8:
                            _id = False

                if not _id and evalue >= minevalue:
                    # get reading frame
                    hits[strand[name]] += 1
                    # poscount[rf] = poscount[rf].union(set(range(start,stop)))
                    for i in range(start, stop):
                        if i not in posgenes[name][strand[name]]:
                            posgenes[name][strand[name]][i] = 1
                        else:
                            posgenes[name][strand[name]][i] += 1

            # geting best rf
            strandmax = -1

            for strandx in hits:
                if hits[strandx] > hits[strandmax]:
                    strandmax = strandx
            strand[name] = strandmax
            posgenes[name] = posgenes[name][strandmax]

            todel = []
            # Cut the hits
            if percent:
                endcut = hits[strandmax] * cutoff
                for i in posgenes[name]:
                    if posgenes[name][i] < endcut:
                        todel.append(i)
            else:
                for i in posgenes[name]:
                    if posgenes[name][i] < cutoff:
                        todel.append(i)
            for i in todel:
                del posgenes[name][i]
            # geting list of the Hits
            posgenes[name] = [i for i in posgenes[name]]
            posgenes[name].sort()

        except Exception:
            logging.warning("blastfile error: could not open file %s\n" % _file)

    featurelist = []
    # For every Gene
    for name in posgenes:
        if len(posgenes[name]) > 0:
            hits = []
            start = posgenes[name][0]
            lasti = posgenes[name][0] - 1
            # get singele hits for one gene
            for i in posgenes[name]:
                if lasti != i - 1:
                    hits.append([start, lasti])
                    start = i
                lasti = i
            hits.append([start, lasti])

            for hit in hits:
                featurelist.append(
                    feature(
                        name,
                        type_from_name(name),
                        int(hit[0]),
                        int(hit[1]),
                        int(strand[name]),
                        "blastn",
                    )
                )

    return featurelist


def _clip_predictions(h1, h2, idx, values):
    """
    function for clipping two overlapping predictions
    intended use: predictions with same name and reasonable height difference

    @param[in] h1 a prediction
    @param[in] h2 another prediction
    @param[in] idx position index for the following three arrays idx[name][rf][p] gives the array position for hit at p in rf
    @param[in] posgenes the number of hits at a position
    @param[in] value    the average evalue at a position
    @param[in] query    the relative query position at a position
    @param circular circularity
    @param length length of the mitogenome

    @return if one hit is inlcuded in the other then then (None, hit2),
        otherwise (hit1,hit2). Thus it has to be taken care that the prefered hit is hit1!
    """
    swap = False

    #     name = h1.name

    # h1 starts left of h2
    if h1.start < h2.start:
        # h1 ---------
        # h2   ----
        if h2.stop <= h1.stop:
            return None, h2  # return h2 in any case, since it is preferred
        # h1 ---------
        # h2      --------
        else:
            swap = False
            hit1, hit2 = h1, h2

    # h2 starts left of h1
    else:
        # h1   ----
        # h2 ---------
        if h1.stop <= h2.stop:
            return None, h2  # return h2 in any case, since it is preferred
        # h1      --------
        # h2 ---------
        else:
            swap = True
            hit1, hit2 = h2, h1

    #     if h1.start <= h2.start <= h1.stop <= h2.stop:
    #         swap = False
    #         hit1, hit2 = h1, h2
    #     elif h2.start <= h1.start <= h2.stop <= h1.stop:
    #         swap = True
    #         hit1, hit2 = h2, h1
    #     else:
    # #        sys.stderr.write( "clip predictions: got strange predictions\n %s %s" % ( str( h1 ), str( h2 ) ) )
    #         return None, h2

    logging.debug("clip %s %s" % (hit1, hit2))

    # get step size
    if abs(hit1.rf) == 4:
        stp1 = 1
    else:
        stp1 = 3

    if abs(hit2.rf) == 4:
        stp2 = 1
    else:
        stp2 = 3

    # move hit2.start to the right as long as the positions are worse than those of hit1
    # move hit1.start to the left as long as the positions are worse than
    # those of hit2
    for i in range(hit2.start, min(hit2.stop, hit1.stop + 1), stp2):
        if values[idx[hit2.rf][i]].get_score() <= values[idx[hit1.rf][i]].get_score():
            hit2.start = i
        else:
            break

    for j in range(hit1.stop, max(hit1.start, hit2.start - 1), -1 * stp1):
        if values[idx[hit1.rf][j]].get_score() <= values[idx[hit2.rf][j]].get_score():
            hit1.stop = j
        else:
            break

    # if there is still overlap -> move the endposition of the feature that is
    # worse in sum
    # var to take care to equally remove from features of the same height
    i = 0
    while (
        hit1.start <= hit1.stop
        and hit2.start <= hit2.stop
        and feature.cap(hit1.start, hit1.stop, hit2.start, hit2.stop, False, 0) > 0
    ):
        #         print "tt", hit1, hit2
        if (i % 2) == 0 and hit1.score >= hit2.score:
            hit1.stop = hit1.stop - stp1

        if (i % 2) == 1 and hit2.score >= hit1.score:
            hit2.start = hit2.start + stp2

        i += 1

    logging.debug("clip %s %s" % (hit1, hit2))

    # recalculate the hits, i.e. get new position sets and re-initialize
    if hit1.start > hit1.stop:
        hit1 = None
    else:
        try:
            _update_hit(hit1, idx, values)
        except UpdateZeroScore:
            hit1 = None

    if hit2.start > hit2.stop:
        hit2 = None
    else:
        try:
            _update_hit(hit2, idx, values)
        except UpdateZeroScore:
            hit2 = None

    #     positions = []
    #     for i in xrange( hit1.start, hit1.stop + 1, 1 ):
    #         positions.append( values[idx[name][name][hit1.rf][i]] )
    #     if positions == []:
    #         hit1 = None
    #     else:
    #         hit1 = feature.blast_feature( hit1.name, hit1.type, hit1.strand, positions, mito = 1, rf = hit1.rf )
    #
    #     positions = []
    #     for i in xrange( hit2.start, hit2.stop + 1, 1 ):
    #         positions.append( values[idx[name][hit2.rf][i]] )
    #     if positions == []:
    #         hit2 = None
    #     else:
    #         hit2 = feature.blast_feature( hit2.name, hit2.type, hit2.strand, positions, mito = 1, rf = hit2.rf )

    if hit1 is None:
        pass
    elif hit2 is None:
        hit1, hit2 = hit2, hit1
    elif hit1 is None and hit2 is None:
        logging.error("clip predictions yielded 2 None features")
        raise Exception("clip predictions yielded 2 None features")
    elif swap:
        hit1, hit2 = hit2, hit1

    return hit1, hit2


def _create_initial_predictions(idx, values, cutoff):
    """
    Determine features, i.e. consecutive stretches of positions (in the
    same reading frame) which have blast hits. In this version a maximum
    (and threshold) is computed and the threshold is applied per stretch.

    This is done separately per gene and reading frame.
    For every hit the average (over the length of the feature)
    evalue and height (number of hits) is computed.
    Furthermore the average query start and stop position is stored.query,
    @param idx: idx[name][rf][pos] gives the idx in the values array where the data for
        gene "name" in reading frame "rf" at position "pos" can be found
    @param values a list of position_values (which contain the position and summed scores, evalues etc)
    @param cutoff cutoff cutoff percentage

    @return: a list of (blast_)features
    """

    hits = []  # the list of features for each gene

    # do for each gene and reading frame:
    for name in idx:
        for rf in idx[name]:
            conspos = []  # sets of consecutive positions
            conspos.append([])
            positionlist = [p for p in list(idx[name][rf].keys())]
            positionlist.sort()
            if rf > 0:
                strand = 1
            else:
                strand = -1

            # -1 -> then it works for the 1st feature
            lastpos = positionlist[0] - 1
            for i in positionlist:
                # if continuation of last position -> do nothing
                if lastpos == i - 1:
                    pass
                # if disrupted:
                # -> finish old feature
                # -> start a new one
                else:
                    conspos.append([])
                p = idx[name][rf][i]
                conspos[-1].append(values[p])
                lastpos = i

            for cp in conspos:
                # determine maximum position

                maxp = max(cp, key=lambda x: x.get_score())
                maxv = maxp.get_score() * cutoff
                maxi = cp.index(maxp)

                #                 logging.warn( "%s [%d,%d]" % ( name, cp[0].get_position(), cp[-1].get_position() ) )
                #                 logging.warn( " ".join( [str( x.get_position() ) for x in cp] ) )
                #                 logging.warn( "maxp %s" % str( maxp ) )
                #                 logging.warn( "maxv %s %s" % ( str( maxp.get_score( ) ), str( maxv ) ) )
                #                 logging.warn( "maxi %s" % str( maxi ) )

                left = maxi
                right = maxi
                while left > 0:
                    #                     logging.warn( "check %d %f" % ( l - 1, cp[l - 1].get_score( ) ) )
                    if cp[left - 1].get_score() >= maxv:
                        left -= 1
                    else:
                        break
                while right < (len(cp) - 1):
                    #                     logging.warn( "check %d %f" % ( r + 1, cp[r + 1].get_score( ) ) )
                    if cp[right + 1].get_score() >= maxv:
                        right += 1
                    else:
                        break
                #                 logging.warn( "===> %d-%d" % ( l, r ) )
                #                 logging.warn( " ".join( [str( x.get_position() ) for x in cp[l:r + 1]] ) )
                hits.append(
                    feature.blast_feature(
                        name,
                        type_from_name(name),
                        strand,
                        cp[left : right + 1],
                        mito=1,
                        rf=rf,
                    )
                )
    return hits


def _create_initial_predictions_globth(idx, values, cutoff):
    """
    this is the old version applying the threshold per reading frame

    determine features, i.e. consecutive stretches of positions (in the
    same reading frame) which have blast hits
    this is done separately per gene and reading frame
    for every hit the average (over the length of the feature)
    evalue and height (number of hits) is computed
    furthermore the average query start and stop position is stored.query,
    @param idx: idx[name][rf][pos] gives the idx in the values array where the data for
        gene "name" in reading frame "rf" at position "pos" can be found
    @param values a list of position_values (which contain the position and summed scores, evalues etc)
    @param cutoff cutoff cutoff percentage

    @return: a list of (blast_)features
    """

    hits = []  # the list of features for each gene

    conspos = []  # a list of consecutive positions

    # do for each gene and reading frame:
    for name in idx:
        for rf in idx[name]:
            # get the maximum quality
            mx = 0
            for p in idx[name][rf]:
                i = idx[name][rf][p]
                mx = max(mx, values[i].get_score())
            endcut = mx * cutoff

            # get all positions that have a quality >= x * maximum quality
            positionlist = [
                p
                for p in list(idx[name][rf].keys())
                if values[idx[name][rf][p]].get_score() >= endcut
            ]
            positionlist.sort()

            if rf > 0:
                strand = 1
            else:
                strand = -1

            p = idx[name][rf][positionlist[0]]
            conspos.append(values[p])
            # -1 -> then it works for the 1st feature
            lastpos = positionlist[0] - 1
            #            nfeat = feature.blast_feature( name, type_from_name(name), \
            #                        positionlist[0], positionlist[0] - 1, \
            #                        strand, "mitos", rf = rf,
            #                        qstart = values[p].query, qstop = 0, score = 0 )

            for i in positionlist:
                # if continuation of last position -> do nothing
                if lastpos == i - 1:
                    pass
                # if disrupted:
                # -> finish old feature
                # -> start a new one
                else:
                    #                    if nfeat.start != nfeat.stop:
                    # #                        nfeat.evalue /= ( nfeat.stop - nfeat.start )
                    # #                        nfeat.height /= float( nfeat.stop - nfeat.start )
                    #
                    #                        pstart = idx[name][rf][nfeat.start]
                    #                        pstop = idx[name][rf][nfeat.stop]
                    #                        nfeat.qstop = values[pstop].query
                    #                        if values[pstart].query > values[pstop].query:
                    #                            nfeat.qstart, nfeat.qstop = nfeat.qstop, nfeat.qstart
                    #                        hits.append( nfeat )
                    hits.append(
                        feature.blast_feature(
                            name, type_from_name(name), strand, conspos, mito=1, rf=rf
                        )
                    )
                    conspos = []

                #                    nfeat = feature.blast_feature( name, type_from_name(name), i, i, \
                #                                strand, "mitos", rf = rf, evalue = 0,
                # qstart = values[p].query, qstop = 0, height = 0 )

                p = idx[name][rf][i]
                # continue counting
                conspos.append(values[p])
                #                nfeat.bitscore += values[p].bitscore
                #                nfeat.evalue += values[p].evalue
                #                nfeat.height += values[p].height
                lastpos = i

            # finish last feature
            hits.append(
                feature.blast_feature(
                    name, type_from_name(name), strand, conspos, mito=1, rf=rf
                )
            )
            conspos = []
    #                nfeat.evalue /= ( nfeat.length( False, 0 ) )
    #                nfeat.height /= float( nfeat.length( False, 0 ) )
    #                pstart = idx[name][rf][nfeat.start]
    #                pstop = idx[name][rf][nfeat.stop]
    #                nfeat.qstop = values[pstop].query
    #                if values[pstart].query > values[pstop].query:
    #                    nfeat.qstart, nfeat.qstop = nfeat.qstop, nfeat.qstart
    #                hits.append( nfeat )
    return hits


def findstart(posi, code, sequence, length, strand, aarng):
    """
    search a position in the close neighbourhood of posi
    that is a start position

    @param posi current start position
    @param code genetic code id
    @param sequence sequence (doubled for circular=True)
    @param length the length of the sequence (undoubled)
    @param strand strand
    @param length length of the sequence
    @param aarng how many aminoacids before and after a start should be searched
    @return a better position or the original one
    """

    #    start, stop = codon.parsecodon( code )
    start = CodonTable.unambiguous_dna_by_id[code].start_codons
    # stop = CodonTable.unambiguous_dna_by_id[code].stop_codons
    #     logging.debug( "find_start %d %d" % ( posi, strand ) )

    # search from the inside of the gene to the outside for a stop codon
    # if one is found then start the search for a start codon behind this
    # stop codon
    if strand == 1:
        frm = posi - aarng * 3
        rng = list(
            range(min(length - 3, posi + aarng * 3), max(0, posi - aarng * 3 - 1), -3)
        )
    else:
        frm = posi + aarng * 3
        rng = list(
            range(max(2, posi - aarng * 3), min(length - 1, posi + aarng * 3 + 1), 3)
        )

    #     logging.debug( "search stop @ start rng %s" % ( str( rng ) ) )
    #     for i in rng:
    # #         logging.debug( "? stop %d" % ( i ) )
    #         if strand == 1:
    #             cdn = sequence.subseq( i, i + 2, strand )
    #         else:
    #             cdn = sequence.subseq( i - 2, i, strand )
    #
    #         if str(cdn) in stop:
    #             if strand == 1:
    #                 frm = i + 3
    #             else:
    #                 frm = i - 3
    #             logging.debug( "found stop %s -> frm %d" % ( cdn, frm ) )
    #             break

    # search for a start codon
    # first the outside is scanned then the inside
    # each of the two scans starts from positions close to the original
    # position
    if strand == 1:
        rng = list(range(max(0, frm), min(length - 3, posi + aarng * 3 + 1), 3))
        rng.sort(key=lambda x: (x > posi, abs(x - posi)))
    else:
        rng = list(range(min(length - 1, frm), max(2, posi - aarng * 3 - 1), -3))
        rng.sort(key=lambda x: (x < posi, abs(posi - x)))

    logging.debug("search start from %d rng %s" % (frm, str(rng)))
    for i in rng:
        #         logging.debug( "? start %d" % ( i ) )
        if strand == 1:
            cdn = sequence.subseq(i, i + 2, strand)
        else:
            cdn = sequence.subseq(i - 2, i, strand)
        logging.debug("%d %s" % (i, cdn))

        if str(cdn) in start:
            logging.debug("found start %s -> %d" % (cdn, i % length))
            return i % length

    if strand == 1:
        return (max(frm, posi)) % length
    else:
        return (min(frm, posi)) % length


def findstop(posi, code, sequence, length, strand, aarng, maxdegenerated=0):
    """
    search a position in the close neighbourhood of posi
    that is a stop position

    @param posi current stop position
    @param code genetic code id
    @param sequence sequence (doubled for circular=True)
    @param length the length of the sequence (undoubled)
    @param strand strand
    @param aarng how many aminoacids before and after a start should be searched
    @param maxdegenerated maximum number of degenerated positions allowed [0,1,2]

    @return a better position or the original one
    """
    if 0 > maxdegenerated or maxdegenerated > 2:
        Exception("Invalid stop codon degeneration %s" % maxdegenerated)

    for degenerated in range(0, maxdegenerated + 1):
        # get stop codons and degenerate them
        # XYA -> XYN
        # XAA -> XNN
        stop = list(CodonTable.unambiguous_dna_by_id[code].stop_codons)
        for j in range(len(stop)):
            st = [x for x in stop[j]]
            for k in range(len(st) - 1, len(st) - degenerated - 1, -1):
                if st[k] != "A":
                    break
                st[k] = "N"
            stop[j] = trna.codon("".join(st), "codon")

        if strand == 1:
            rng = list(
                range(
                    max(0, posi - aarng * 3), min(length - 3, posi + aarng * 3 + 1), 3
                )
            )
        else:
            rng = list(
                range(
                    min(length - 1, posi + aarng * 3), max(2, posi - aarng * 3 - 1), -3
                )
            )

        logging.debug("search stop rng %s" % (str(rng)))
        for i in rng:
            if (strand == 1 and i < 2) or (strand == -1 and i >= length - 2):
                continue

            if strand == 1:
                cdn = sequence.subseq(i - 2, i, strand)
            else:
                cdn = sequence.subseq(i, i + 2, strand)

            # check if codon is a stop codon
            for s in stop:
                if cdn.isequal(s):
                    logging.debug(
                        "found stop %s -> %d (%d)"
                        % (cdn, i, (i - (strand * degenerated)) % length)
                    )
                    # return position i
                    # - degenerated positions
                    # - modulo sequence length -> i\in [0:length] note this holds
                    #   automatically for linear sequences
                    if degenerated > 0:
                        logging.error("found degenerated stop %s" % (cdn))
                    return (i - (strand * degenerated)) % length

    return posi


def _aggregate_iter(blasthits, values, featlen, gname=""):
    """
    the actual function to get the aggregation of the blast hits

    this is the slow [but simple] implementation iterating for each
    blast hit over the included positions

    is still included as reference implementation

    @param[in] blasthits a (nonempty) list of blast hits to aggregate
    @param[in] length the length of the (original) sequence (only important for circular seq)
    @param[out] values the aggregated values (evalue, bitscore, avg. relative query position)
    @param[in] featlen dictionary containing the lengths of the query sequences
    @param[in] gname gene name (only needed for debugging)
    @return a mapping from genomic positions to the index in values
        i.e. idx[rf][i] gives the index for genomic position i in reading frame rf
    """

    idx = {
        1: dict(),
        2: dict(),
        3: dict(),
        4: dict(),
        -1: dict(),
        -2: dict(),
        -3: dict(),
        -4: dict(),
    }
    blasthits.sort(key=lambda x: (x.start, x.stop))
    # this is the simple (but slow) reference implementation of the hit agglomeration.
    # - for a blast hit for each position in the range (start .. end) covered by the hit
    #   the values (evalue, bitscore, relative query position) are summed up
    # - this is done for each blast hit.
    #     cnt = 0
    for h in blasthits:
        #         logging.debug( "H %s" % str( h ) )
        flen = featlen[gname][h.qname]

        #         fac = fractions.Fraction(
        #                 feature.length( min( b.qstart, b.qstop ), max( b.qstart, b.qstop ), False, 0 ) , \
        #                 b.length( circular, length ) )
        #         fac = decimal.Decimal( feature.length( min( h.qstart, h.qstop ), max( h.qstart, h.qstop ), False, 0 ) ) / decimal.Decimal( h.length( False, 0 ) )
        fac = feature.length(
            min(h.qstart, h.qstop), max(h.qstart, h.qstop), False, 0
        ) / float(h.length(False, 0))
        for i in feature.crange(h.start, h.stop + 1, 1, False, 0):
            if i not in idx[h.rf]:
                idx[h.rf][i] = len(values)
                values.append(position_values(i))
            p = idx[h.rf][i]
            # - why length-1: the 1st position is at distance 0
            # - all positions computed wrt start (in the reference)
            if h.strand == 1:
                relquery = (
                    h.qstart + fac * (feature.length(h.start, i, False, 0) - 1)
                ) / float(flen)
                relyreuq = (
                    (flen - h.qstart) - fac * (feature.length(h.start, i, False, 0) - 1)
                ) / float(flen)
            else:
                relquery = (
                    h.qstop - fac * (feature.length(h.start, i, False, 0) - 1)
                ) / float(flen)
                relyreuq = (
                    (flen - h.qstop) + fac * (feature.length(h.start, i, False, 0) - 1)
                ) / float(flen)
            values[p].add(h.evalue, h.bitscore, relquery, relyreuq)
    #         if cnt > 1:
    #             break
    #         else:
    #             cnt += 1

    return idx


def _aggregate_sweep(blasthits, values, featlen, gname):
    """
    the actual function to get the aggregation of the blast hits

    this is the fast (line sweep) implementation

    @param[in] blasthits a (nonempty) list of blast hits to aggregate
    @param[in] circular true: the sequence is treated as circular
    @param[in] length the length of the (original) sequence (only important for circular seq)
    @param[out] values the aggregated, i.e. summed, values (evalue, bitscore,
        relative query position)
    @param[in] featlen dictionary containing the lengths of the query sequences
    @param[in] gname gene name
    @return a mapping from genomic positions to the index in values
        i.e. idx[rf][i] gives the index for genomic position i in reading frame rf
    """

    idx = {
        1: dict(),
        2: dict(),
        3: dict(),
        4: dict(),
        -1: dict(),
        -2: dict(),
        -3: dict(),
        -4: dict(),
    }

    # the idea of this function is to sweep once over the genome to compute
    # the values (quality, rel.query position ...) for each position
    #
    # the key idea is the observation that the aggregated values either
    # change only at the begin or end of a blast hit (height, evalue,
    # bitscore) or additionally at each position by a constant value
    # (relative query position)
    #
    # so we get a list of events (the set of start and end positons)
    # which point to the corresponding blast hits that start or end there
    #
    # at each event
    # - the values from the last event to the current event are actually
    #   added to the values array.
    #   the relative query position depends on the distance from the start
    #   or stop. therefore a constant value is added for each position to
    #   the current values
    # - the aggregated values are updated (add/remove contribution
    #   of hits starting/ending here)

    #     for b in blasthits:
    #         if b.strand == 1:
    #             logging.debug( "%s %d %d" % ( b, b.qstart, b.qstop ) )
    #         else:
    #             logging.debug( "%s %d %d" % ( b, b.qstop, b.qstart ) )

    # initialise the list of events
    blasthits.sort(key=lambda x: (x.start, x.stop))
    events = {}
    #     cnt = 0
    for b in blasthits:
        #         logging.debug( "H %s" % str( b ) )
        try:
            events[b.start].append(b)
        except KeyError:
            events[b.start] = [b]

        try:
            events[b.stop + 1].append(b)
        except KeyError:
            events[b.stop + 1] = [b]
    #         if cnt > 1:
    #             break
    #         else:
    #             cnt += 1
    sortedevents = sorted(events.keys())
    # initialize the current values dict (cv) these are the aggregated
    # values that are updated at each event. in addition the constant
    # values (dq) to be added at each position to the relative queries
    # are initialised. both values are stored separately per reading
    # frame. if there are features spaning 0 then there is a virtual event
    # added (by setting prevev to 0)

    po = blasthits[0].start  # current position
    cv = {
        1: position_values(po),
        2: position_values(po),
        3: position_values(po),
        4: position_values(po),
        -1: position_values(po),
        -2: position_values(po),
        -3: position_values(po),
        -4: position_values(po),
    }
    dq = {1: 0, 2: 0, 3: 0, 4: 0, -1: 0, -2: 0, -3: 0, -4: 0}

    # processing will start at the first event (the one with smallest position)
    prevev = None

    #     logging.debug( "%s init %s" % ( gname, str( cv ) ) )
    for e in sortedevents:
        #         logging.debug( "%s event %s pevent %s" % ( gname, str( e ), str( prevev ) ) )
        #         logging.debug( "%s    cv %s" % ( gname, str( cv ) ) )
        #         logging.debug( "%s    dq %s" % ( gname, str( dq ) ) )

        # add sum for everything from the last event to the current
        # if there is something to add (height > 0)
        # and update the aggregated relative query position by dq
        if prevev is not None:
            for f in cv:
                if cv[f].get_height() <= 0:
                    continue
                for p in feature.crange(prevev, e, 1, False, 0):
                    cv[f].set_position(p)
                    idx[f][p] = len(values)
                    values.append(copy.deepcopy(cv[f]))
                    #                     logging.debug( "%s add %d %s " % ( gname, f, str( values[-1] ) ) )
                    #                     cv[f].set_query( ( cv[f].get_query() ) )
                    #                     cv[f].set_yreuq( ( cv[f].get_yreuq() ) )
                    #                     if p != prevev:
                    cv[f].set_query((cv[f].get_query() + dq[f]))
                    cv[f].set_yreuq((cv[f].get_yreuq() - dq[f]))

        prevev = e
        # add contribution of hits starting here
        for b in events[e]:
            if e != b.start:
                continue
            flen = featlen[gname][b.qname]
            fac = (
                feature.length(min(b.qstart, b.qstop), max(b.qstart, b.qstop), False, 0)
                / float(b.length(False, 0))
            ) / float(flen)
            dq[b.rf] += b.strand * fac
            if b.strand == 1:
                relquery = b.qstart / float(flen)
                relyreuq = (flen - b.qstart) / float(flen)
            else:
                relquery = (b.qstop) / float(flen)
                relyreuq = (flen - b.qstop) / float(flen)
            cv[b.rf].add(b.evalue, b.bitscore, relquery, relyreuq)

        # remove contribution of hits ending here
        for b in events[e]:
            if e == b.start:
                continue
            flen = featlen[gname][b.qname]
            fac = feature.length(
                min(b.qstart, b.qstop), max(b.qstart, b.qstop), False, 0
            ) / float(b.length(False, 0))
            dq[b.rf] -= (b.strand * fac) / float(flen)
            if b.strand == 1:
                #                 relquery = b.qstart / float( flen )
                #                 relyreuq = ( flen - b.qstart ) / float( flen )
                relquery = (
                    b.qstart + fac * (feature.length(b.start, b.stop, False, 0))
                ) / float(flen)
                relyreuq = (
                    (flen - b.qstart)
                    - fac * (feature.length(b.start, b.stop, False, 0))
                ) / float(flen)
            else:
                # relquery = ( b.qstop ) / float( flen )
                # relyreuq = ( flen - b.qstop ) / float( flen )
                relquery = (
                    b.qstop - fac * (feature.length(b.start, b.stop, False, 0))
                ) / float(flen)
                relyreuq = (
                    (flen - b.qstop) + fac * (feature.length(b.start, b.stop, False, 0))
                ) / float(flen)
            cv[b.rf].sub(b.evalue, b.bitscore, relquery, relyreuq)

    return idx


def _blasthandle(
    blastfile: str,
    minevalue: float,
    values: List[position_values],
    featlen: Dict[str, Dict[str, int]],
    prot: bool = True,
    ftype: Literal["xml", "tab"] = "tab",
):
    """
    read a blast file

    for each position the following values are determined:
    - the number of blast hist including the position
    - the average e-value of these blast hits and
    - the average relative query position of these hits

    @param[in] blastfile a file name of a tabular (-m 8) or xml  (-m 7) blast output
    @param[in] minevalue e-value threshold to be applied, i.e. smaller evalues
        are discarded
    @param[out] values aggregated values for each position (number of hits for
        the position, average of the e-value and bitscore, and average of the
        relative query positions)
    @param[in] featlen dictionary containing the lengths of the query sequences
    @paramp[in] prot
        true: the blast results are for proteins -> will be assigned to reading frames [-3,-2,-1,1,2,3]
        false: nucleotide matches -> reading frames 4, -4 depending on the strand
    @param[in] ftype filetype tab: tabular BLAST output, xml: XML output

    @return a mapping from genomic positions to the index in values
    i.e. idx[rf][i] gives the index for genomic position i in reading frame rf
    """

    gname = blastfile.split(".")[-2]  # gene name

    # iterate through the blast file lines (blast hits) for the current gene
    blasthandle = open(blastfile, "r")
    #     blasthits = {}  # separated list of blast hits per query sequence
    blasthits = []
    if ftype == "xml":
        handle = Bio.Blast.NCBIXML.parse(blasthandle)
        for fh in handle:
            for alignmen in fh.alignments:
                qname = alignmen.hit_def.split()[0]

                for hsp in alignmen.hsps:
                    if 0 == hsp.expect:
                        evalue = 100.0
                    else:
                        evalue = -1.0 * math.log10(hsp.expect)

                    if evalue < minevalue:
                        continue

                    start = int(hsp.query_start) - 1
                    stop = int(hsp.query_end) - 1

                    querystart = int(hsp.sbjct_start) - 1
                    querystop = int(hsp.sbjct_end) - 1
                    if prot:
                        querystart *= 3
                        querystop *= 3

                    bitscore = hsp.bits

                    if prot:
                        rf = hsp.frame[0]
                        if start > stop:
                            raise Exception("start > stop in blastx results")
                    else:
                        if start > stop and querystart < querystop:
                            rf = -4
                        elif start < stop and querystart > querystop:
                            rf = -4
                        else:
                            rf = 4

                    if start > stop:
                        start, stop = stop, start
                    if querystart > querystop:
                        querystart, querystop = querystop, querystart

                    if rf < 0:
                        strand = -1
                    else:
                        strand = 1

                    #                     if not alignmen.accession in blasthits:
                    #                         blasthits[alignmen.accession] = []

                    # add feature (make end position excluded [they are alreadt
                    # 0 based])
                    nfeat = feature.blasthit_feature(
                        gname,
                        type_from_name(gname),
                        qname,
                        start,
                        stop + 1,
                        strand,
                        rf,
                        querystart,
                        querystop,
                        evalue,
                        bitscore,
                    )
                    #                     blasthits[alignmen.accession].append( nfeat )
                    blasthits.append(nfeat)

    elif ftype == "tab":
        for line in blasthandle:
            # split the line in the single informations
            cols = line.split()

            if not cols or len(cols) < 12:
                continue

            qname = cols[-11]  # query name

            # getting evalue for the cut of evalue
            if cols[-2] == "0.0":
                evalue = 200.0
            else:
                evalue = -1.0 * math.log10(float(cols[-2]))

            if evalue <= 0.0 or evalue < minevalue:
                continue

            start = int(cols[-6]) - 1
            stop = int(cols[-5]) - 1
            querystart = int(cols[-4]) - 1
            querystop = int(cols[-3]) - 1
            if prot:
                querystart *= 3
                querystop *= 3

            bitscore = float(cols[-1])

            # get strand
            if start > stop and querystart < querystop:
                strand = -1
            elif start < stop and querystart > querystop:
                strand = -1
            elif start < stop and querystart < querystop:
                strand = 1
            elif start > stop and querystart > querystop:
                strand = 1
            else:
                logging.error("[Blast] error: start = stop")
                sys.exit()

            if start > stop:
                start, stop = stop, start
            if querystart > querystop:
                querystart, querystop = querystop, querystart
            # get reading frame
            if prot:
                rf = (strand) * (start % 3 + 1)
            else:
                rf = 4 * strand

            #             if not query in blasthits:
            #                 blasthits[query] = []

            # add feature (make end position excluded [they are alreadt 0
            # based])
            nfeat = feature.blasthit_feature(
                gname,
                type_from_name(gname),
                qname,
                start,
                stop,
                strand,
                rf,
                querystart,
                querystop,
                evalue,
                bitscore,
            )
            blasthits.append(nfeat)
    #             blasthits[query].append( nfeat )
    else:
        raise Exception("Unknown BLAST file type")

    blasthandle.close()

    if len(blasthits) == 0:
        return dict()

    # aggregate blast hits
    #     idx = _aggregate_iter( blasthits, values, featlen, gname )
    idx = _aggregate_sweep(blasthits, values, featlen, gname)

    # compute average values of query positions and e-values
    # and delete empty reading frames
    for rf in list(idx.keys()):
        #         for i in sorted( idx[rf].iterkeys() ):
        #             logging.debug( "%s %d %s " % ( gname, rf, values[idx[rf][i]] ) )
        # #             values[idx[rf][i]].finalize()
        for p in list(idx[rf].values()):
            values[p].finalize()
        if len(idx[rf]) == 0:
            del idx[rf]

    return idx


def _get_codon(p, strand, sequence):
    """
    code to obtain the codon starting/ending at a postion
    in the sequence on the +/- strand
    @param[in] p the position
    @param[in] strand the strandedness
    @param[in] circular treat sequence as circular
    @param[in] sequence the sequence
    @return the codon as codon object
    """
    cdn = sequence.subseq(p, p + 2, strand)
    cdn = str(cdn)
    return trna.codon(cdn, "codon")


def _improve_start_stop(
    hit,
    seq,
    code,
    startcod,
    stopcod,
    len_pval,
    idx,
    values,
    cpyprt,
    forbidden,
    circular,
):
    """
    @param[in,out] the hit to improve, given as list of parts
    @param[in] seq the sequence
    @param[in] startcod start codon stats
    @param[in] stopcod stop codon stats
    @param[in] featlenstats a dictionary containing for each gene the
        mean (mu), standard deviation (sigma), maximum (max), number of data points (n) and pvalue list (pv)
    @param[in] idx index from (rf,position) to index in values
    @param[in] values quality values
    @param[in] cpyprt index of (name,copy) to last part
    @param[in] forbidden a list of features that are not allowed to be overlapped
    @param[in] circular genome circularity
    @return plot data
    """
    if len(hit) > 1:
        logging.debug(
            "_improve_start_stop %s .. %s %s"
            % (str(hit[0]), str(hit[-1]), str(circular))
        )
    else:
        logging.debug("_improve_start_stop %s %s" % (str(hit[0]), str(circular)))
    ret = []

    # sort the parts of the hit by query start position
    hit.sort(key=lambda k: k.qstart)
    if hit[0].strand < 0:
        hit.reverse()

    # save old start and stop position to determine if the
    # hit needs to be updated
    ostart = hit[0].start
    ostop = hit[-1].stop

    # determine maximum and minimum quality mx mn over all parts
    mx = 0
    mn = float("inf")
    for i in range(len(hit)):
        mx = max(mx, max([x.get_score() for x in hit[i]._positions]))
        mn = min(mn, min([x.get_score() for x in hit[i]._positions]))
    #     logging.debug( "MIN %f MAX %f" % ( mn, mx ) )

    if mn < 0:
        raise Exception("_improve_start_stop found negative quality value(s)")

    # find maximum quality plateau plat=(p1,p2) where p1 and p2 mark the boundaries
    # it needs to be consistent with the reading frame
    plat = [hit[0].start, hit[-1].stop - 2]
    rng = [hit[0].start, hit[-1].stop - 2]

    #     for i in hit[0]._positions:
    #         if ( feature.length( ostart, i.get_position(), circular, len( seq ) ) - 1 ) % 3 != 0:
    #             continue
    #         if i.get_score() >= mx:
    #             plat[0] = i.get_position()
    #             break
    #     for i in hit[-1]._positions[::-1]:
    #         if feature.length( i.get_position(), ostop, circular, len( seq ) ) % 3 != 0:
    #             continue
    #         if i.get_score() >= mx:
    #             plat[1] = i.get_position()
    #             break
    # shrink plateau to the middle
    mod = True
    while mod:
        mod = False
        if plat[0] + 3 < hit[0].stop and plat[0] + 3 <= plat[1] - 3:
            plat[0] += 3
            mod = True
        if plat[1] - 3 > hit[-1].start and plat[0] + 3 <= plat[1] - 3:
            plat[1] -= 3
            mod = True
    # logging.debug("\tplat %s %s with value %f " % (plat[0], plat[1], mx))

    #     if plat[1] < plat[0] and not circular:
    #         logging.error( "_improve_start_stop: invalid plat %d %d" % ( plat[0], plat[1] ) )

    # determine stops/forbidden regions preceeding and following the plateu
    # this determines the range to search for the start and stop
    i = plat[0]
    while (not circular and i >= 0) or (
        circular and (i <= plat[0] or (i > plat[0] and i > plat[1]))
    ):
        cdn = _get_codon(i, hit[0].strand, seq)
        if cdn.isstop(code):
            rng[0] = i
            # logging.debug("left stop %d %d due to stop %s" % (i, rng[0], str(_get_codon(i, hit[0].strand, seq))))
            break
        # elif str(cdn) == "NNN":
        #     if i == plat[0]:
        #         rng[0] = plat[0]
        #     else:
        #         rng[0] = (i + 3) % len(seq)
        #     logging.debug("left stop %d %d due to ambig %s" % (i, rng[0], str(_get_codon(i, hit[0].strand, seq))))
        #     break
        trnaovl = False
        for t in forbidden:
            cap = feature.cap(i, plat[0], t.start, t.stop, True, len(seq))
            if cap == t.length(True, len(seq)):
                trnaovl = True
                # logging.debug("left stop %d due to forbidden overlap with %s" % (i, str(t)))
                break
        if trnaovl:
            break

        rng[0] = i
        i -= 3
        if circular:
            i = i % len(seq)

    i = plat[1]
    while (not circular and i <= len(seq) - 3) or (
        circular and (i >= plat[1] or (i < plat[1] and i < plat[0]))
    ):
        cdn = _get_codon(i, hit[-1].strand, seq)
        if cdn.isstop(code):
            # logging.debug("right stop %d due to stop %s" % (i, str(cdn)))
            rng[1] = i
            break
        # elif str(cdn) == "NNN":
        #     if i == plat[1]:
        #         rng[1] = plat[1]
        #     else:
        #         rng[1] = (i - 3) % len(seq)
        #     logging.debug("right stop %d %d due to ambig %s" % (i, rng[1], str(cdn)))
        #     break
        trnaovl = False
        for t in forbidden:
            cap = feature.cap(plat[1], i, t.start, t.stop, True, len(seq))
            if cap == t.length(True, len(seq)):
                # logging.debug("right stop %d due to forbidden overlap with %s" % (i, str(t)))
                trnaovl = True
                break
        if trnaovl:
            break
        rng[1] = i
        i += 3
        if circular:
            i = i % len(seq)

    #     rng[0] = plat[0] % 3
    #     rng[1] = len( seq ) - 3
    #     while rng[1] % 3 != rng[0] % 3:
    #         rng[1] -= 1
    #     for i in feature.crange( plat[0], -1, -3, circular, len( seq ) ):
    #         if _get_codon( i, hit[0].strand, seq ).isstop( code ):
    #             rng[0] = i
    #             break
    #         trnaovl = False
    #         for t in forbidden:
    #             cap = feature.cap( i, plat[0], t.start, t.stop, circular, len( seq ) )
    #             if cap == t.length( circular, len( seq ) ):
    #                 trnaovl = True
    #                 break
    #         if trnaovl:
    #             rng[0] = i
    #             break
    #
    #     for i in feature.crange( plat[1], len( seq ) - 2, 3, circular, len( seq ) ):
    #         if _get_codon( i, hit[-1].strand, seq ).isstop( code ):
    #             rng[1] = i
    #             break
    #         trnaovl = False
    #         for t in forbidden:
    # #                 print t
    #             cap = feature.cap( plat[1], i, t.start, t.stop, circular, len( seq ) )
    #             if cap == t.length( circular, len( seq ) ):
    #                 trnaovl = True
    #                 break
    #         if trnaovl:
    #             rng[1] = i
    #             break
    # logging.debug("\trng %s %s" % (str(rng[0]), str(rng[1])))

    # initialize
    # codp: codon probabilities
    # degmal: degeneration malus
    # dend: distance to the corresponding end
    # cod: the codons needed for debugging only
    # quality: qualities
    # quality indicator: (1.0 if >0 and 0.0 if == 0)
    # lenpv: length p-value for pairs of positions
    crit = {}
    for c in ["codp", "degmal", "cod", "dend", "quality", "iquality", "lenpv"]:
        crit[c] = {}

    # get all positions of the hit that are outside of the plateau and inside the range
    # where positions[0] are those between rng[0] and plat[0]
    # and position[1] are those between plat[1] and rng[1]
    # that is positions are sorted outside-to-inside

    positions = [
        list(feature.crange(rng[0], plat[0] + 1, 3, True, len(seq))),
        list(feature.crange(rng[1], plat[1] - 1, -3, True, len(seq))),
    ]
    # logging.error(f"0: {rng[0]}, {plat[0] + 1} {positions[0]}")
    # logging.error(f"0: {rng[1]}, {plat[1] - 1} {positions[1]}")
    # get normalized quality values, more precicely normalized quality+1

    # get codon probabilities and degeneration malus
    #     logging.debug( "%s" % ( str( startcod.get( hit[0].name, {} ) ) ) )
    for i in [0, -1]:
        for p in positions[i]:
            codon = str(_get_codon(p, hit[i].strand, seq))
            codonlst = list(codon)
            codonlst[2] = "N"
            codoni1 = "".join(codonlst)
            codonlst[1] = "N"
            codoni2 = "".join(codonlst)

            cod = codon
            codp = 0.0
            degmal = 1.0

            if (hit[i].strand == 1 and i == 0) or (hit[i].strand == -1 and i == -1):
                codp = startcod.get(hit[0].name, {}).get(codon, 0.0)
            else:
                if not hit[i].name in stopcod:
                    codp = 0.0
                else:
                    if codon in stopcod[hit[i].name]:
                        codp = stopcod[hit[i].name][codon]
                    elif codoni1 in stopcod[hit[i].name]:
                        degmal = 1 / 3.0
                        codp = stopcod[hit[i].name][codoni1]
                        cod = codoni1
                    elif codoni2 in stopcod[hit[i].name]:
                        degmal = 1 / 12.0
                        codp = stopcod[hit[i].name][codoni2]
                        cod = codoni2
                    else:
                        codp = 0.0

            if codp > 0:
                crit["cod"][p] = cod
                crit["codp"][p] = codp
                crit["degmal"][p] = degmal

    #             logging.debug( "%d %d %s %f %f" % ( i, p, cod, codp, degmal ) )

    # from here consider only potential start / stop positions
    for i in [0, 1]:
        tmp = [x for x in positions[i] if x in crit["codp"]]
        positions[i] = tmp

    #     logging.debug( "positions " + str( positions ) )

    # minq = 1 / float( mx )
    #     minq = mn / 10.0
    for i in [0, -1]:
        for p in positions[i]:
            # if there is no max then all should be 1 -> then the position is
            # decided on the other criteria
            if mx == 0:
                crit["quality"][p] = 1
                continue

            try:
                x = idx[hit[i].rf][p]
                crit["quality"][p] = (1 + values[x].get_score()) / float(mx)
            except KeyError:
                crit["quality"][p] = 1 / float(mx)
            if crit["quality"][p] > 0:
                crit["iquality"][p] = 1.0
            else:
                crit["iquality"][p] = 0.0

    #     dquality = {}  # quality difference
    #     for i in [0, -1]:
    #         lastq = 0
    #         for p in positions[i]:
    #             if lastq < crit["quality"][p]:
    #                 lastq = crit["quality"][p]
    #
    #             dquality[p] = crit["quality"][p] - lastq

    #         for j in xrange( len( positions[i] ) ):
    #             dquality[positions[i][j]] = crit["quality"][positions[i][j]] - lastq
    # #             logging.error( "Q %d %d %d %f %f" % ( i, j, positions[i][j], crit["quality"][positions[i][j]], dquality[positions[i][j]] ) )
    #             if j + 1 < len( positions[i] ) and \
    #                     crit["quality"][positions[i][j]] < crit["quality"][positions[i][j + 1]]:
    #                 lastq = crit["quality"][positions[i][j]]

    # the query position normalization values for the part
    # before and after the plateau
    qm = [None, None]
    if hit[0].strand == 1:
        qm[0] = min([x.get_query() for x in hit[0]._positions])
    else:
        qm[0] = min([x.get_yreuq() for x in hit[0]._positions])

    if hit[-1].strand == 1:
        qm[1] = min([x.get_yreuq() for x in hit[-1]._positions])
    else:
        qm[1] = min([x.get_query() for x in hit[-1]._positions])

    #     for x in hit[0]._positions:
    #         logging.debug( "%s q %f y %f %f" % ( str( x ), x.get_query(), x.get_yreuq(), x.get_query() + x.get_yreuq() ) )
    #     logging.debug( "\tqm %f %f %d" % ( qm[0], qm[1], ( plat[0] + plat[1] ) / 2 ) )

    for i in [0, -1]:
        for p in positions[i]:
            v = None
            try:
                x = idx[hit[i].rf][p]
            except KeyError:
                crit["dend"][p] = 1 - qm[i]
            #                 if i == 0:
            #                     crit["dend"][p] = 1 - min( qm[i], abs( p - hit[i].start ) ) / qm[i]
            #                 else:
            #                     crit["dend"][p] = 1 - min( qm[i], abs( p - hit[i].stop ) ) / qm[i]
            #                 logging.debug( "%d %d %f NA" % ( i, p, crit["dend"][p] ) )
            else:
                if (i == 0 and hit[i].strand == 1) or (i == -1 and hit[i].strand == -1):
                    v = values[x].get_query()
                else:
                    v = values[x].get_yreuq()
                crit["dend"][p] = 1 - v
    #                 logging.debug( "dend %d %s %f" % ( p, str( v ), crit["dend"][p] ) )
    #     for i in [0, 1]:
    #         logging.debug( "\tpos[%d]: %s" % ( i, positions[i] ) )

    # get length pvalues for pairs of positions
    for i in positions[0]:
        for j in positions[1]:
            if len(hit) == 1:
                length = feature.length(i, j, True, len(seq))
            else:
                length = feature.length(i, hit[0].stop, True, len(seq))
                length += feature.length(hit[-1].start, j, True, len(seq))
                for k in range(1, len(hit) - 1):
                    length += feature.length(hit[k].start, hit[k].stop, True, len(seq))

            pvn = len_pval.get(hit[0].name, {})
            try:
                pv = pvn[str(length)]
            except KeyError:
                pv = 0

            crit["lenpv"][(i, j)] = pv

    #     # get ranks of the criteria
    #     ranks = {}
    #     for c in crit:
    #         r = {key: rank for rank, key in enumerate( sorted( set( crit[c].values() ), reverse = True ), 1 )}
    #         ranks[c] = {}
    #         for k in crit[c]:
    #             ranks[c][k] = r[crit[c][k]]
    #         logging.debug( "crit[%s] %s" % ( c, crit[c] ) )
    #         logging.debug( "rank[%s] %s" % ( c, ranks[c] ) )

    # get best
    best = 0.0
    for i in positions[0]:
        for j in positions[1]:
            nstart = i + crit["cod"][i].count("N")
            nstop = j + 2 - crit["cod"][j].count("N")

            #             val = 0
            #             cnt = 0
            #             for c in ["codp", "dend", "degmal", "lenpv"]:
            #                 try:
            #                     val += crit[c][i]
            #                     val += crit[c][j]
            #                     cnt += 2
            #                 except:
            #                     val += crit[c][( i, j )]
            #                     cnt += 1
            # #             logging.debug( "%d, %d" % ( val, cnt ) )
            #             val = float( val ) / float( cnt )

            val = 1
            for c in ["codp", "dend", "degmal", "lenpv"]:
                try:
                    val *= crit[c][i]
                    val *= crit[c][j]
                except KeyError:
                    val *= crit[c][(i, j)]

            ret.append(
                "{gene} {strand} {start} {stop} {pva:e} {va:e} {ci} {qi} {di} {fi} {mi} {cj} {qj} {dj} {fj} {mj}\n".format(
                    gene=hit[0].name,
                    strand=hit[0].strand,
                    start=nstart,
                    stop=nstop + 1,
                    va=val,
                    pva=crit["lenpv"][(i, j)],
                    ci=crit["cod"][i],
                    qi=crit["quality"][i],
                    di=crit["dend"][i],
                    fi=crit["codp"][i],
                    mi=crit["degmal"][i],
                    cj=crit["cod"][j],
                    qj=crit["quality"][j],
                    dj=crit["dend"][j],
                    fj=crit["codp"][j],
                    mj=crit["degmal"][j],
                )
            )
            #
            logging.debug(
                "\tstart %d %s q%f d%f c%f d%f "
                % (
                    i,
                    crit["cod"][i],
                    crit["quality"][i],
                    crit["dend"][i],
                    crit["codp"][i],
                    crit["degmal"][i],
                )
            )
            logging.debug(
                "\tstop  %d %s q%f d%f c%f d%f "
                % (
                    j + 2,
                    crit["cod"][j],
                    crit["quality"][j],
                    crit["dend"][j],
                    crit["codp"][j],
                    crit["degmal"][j],
                )
            )
            logging.debug("\t\t len p-len %f => %f " % (crit["lenpv"][(i, j)], val))

            if val > best:
                best = val
                hit[0].start = nstart
                hit[-1].stop = nstop
                # logging.debug("\t\tUPDATING to %d %d" % (nstart, nstop))

    # if the hit got extended to regions where no support from quality
    # values is available then corresponding values (0) need to be inserted
    # at these positions
    logging.debug("\tbest %f" % (best))
    if len(hit) == 1:
        logging.debug("\t==> %s" % (str(hit)))
        if hit[0].start != ostart or hit[0].stop != ostop:
            _update_hit(hit[0], idx, values, circular, len(seq))
    else:
        logging.debug("\t==> %s .. %s" % (str(hit[0]), str(hit[-1])))
        if hit[0].start != ostart:
            _update_hit(hit[0], idx, values, circular, len(seq))
        if hit[-1].stop != ostop:
            _update_hit(hit[-1], idx, values, circular, len(seq))

    #     for i in [0, -1]:
    #         for p in positions[i]:
    #             if crit["codp"][p] == 0:
    #                 continue
    #
    #             if i == 0:
    #                 reg = "start"
    #             else:
    #                 reg = "stop"
    #
    #             ret.append( "{gene} {strand} {genestart} {genestop} {region} {i} {ci} {fi} {qi} {di} {mi}\n".format(
    #                            gene = hit[0].name, strand = hit[0].strand,
    #                            genestart = hit[0].start, genestop = hit[-1].stop,
    #                            region = reg, i = p,
    # ci = crit["cod"][p], fi = crit["codp"][p], qi = quality[p], di =
    # crit["dend"][p], mi = degmal[p] ) )

    for i in range(len(ret)):
        ret[i] = "%d %d %s" % (hit[0].start, hit[-1].stop, ret[i])

    return "".join(ret)


def _improve_start_stop_old(hit, code, sequence, length, idx, values, ststrange=6):
    """
    improve start / stop position (MITOS v1)

    @param hit
    @param sequence
    @param length
    @param ststrange
    @param cpyprt
    """
    logging.debug("adapt stst %s" % (str(hit)))

    if hit[0].strand == 1:
        s = findstart(
            int(hit[0].start), code, sequence, length, hit[0].strand, ststrange
        )
        if s <= hit[0].stop:
            hit[0].start = s
    else:
        s = findstart(
            int(hit[0].stop), code, sequence, length, hit[0].strand, ststrange
        )
        if hit[0].start <= s:
            hit[0].stop = s

    if hit[-1].strand == 1:
        s = findstop(
            int(hit[-1].stop), code, sequence, length, hit[-1].strand, ststrange
        )
        if hit[-1].start <= s:
            hit[-1].stop = s
    else:
        s = findstop(
            int(hit[-1].start), code, sequence, length, hit[-1].strand, ststrange
        )
        if s <= hit[-1].stop:
            hit[-1].start = s

    _update_hit(hit[0], idx, values)

    if len(hit) > 1:
        _update_hit(hit[-1], idx, values)

    logging.debug("=> %s" % (str(hit)))


#     # improve start / stop codons
#     if code is not None:
#         for name in predictions.keys():
#             if name in ["OH", "OL"]:
#                 continue
#             gncstart = startcod[name]
#
#             i = 0
#             while i < len( predictions[name] ):
#                 logging.debug( "adapt stst %s" % ( str( predictions[name][i] ) ) )
#                 while len( predictions[name][i] ) > 0:
# #                     delta = predictions[name][i][0].length( circular, length )
#                     if predictions[name][i][0].strand == 1:
#                         predictions[name][i][0].start = findstart( int( predictions[name][i][0].start ), code, sequence, length, predictions[name][i][0].strand, ststrange, gncstart )
#                     else:
#                         predictions[name][i][0].stop = findstart( int( predictions[name][i][0].stop ), code, sequence, length, predictions[name][i][0].strand, ststrange, gncstart )
# #                     # check if stop overtakes start
#                     if predictions[name][i][0].stop < predictions[name][i][0].start:
#                         logging.debug( "empty %d %d" % ( predictions[name][i][0].start, predictions[name][i][0].stop ) )
#                         del predictions[name][i][0]
#                     else:
#                         break
# #                     # 0 == -1
# #                     if len( predictions[name][i] ) < 1:
# #                         continue
#                 while len( predictions[name][i] ) > 0:
# #                     delta = predictions[name][i][0].length( circular, length )
#                     if predictions[name][i][-1].strand == 1:
#                         predictions[name][i][-1].stop = findstop( int( predictions[name][i][-1].stop ), code, sequence, length, predictions[name][i][-1].strand, ststrange, maxdegenerate )
#                     else:
#                         predictions[name][i][-1].start = findstop( int( predictions[name][i][-1].start ), code, sequence, length, predictions[name][i][-1].strand, ststrange, maxdegenerate )
#                     if predictions[name][i][-1].stop < predictions[name][i][-1].start:
#                         logging.debug( "empty %d %d" % ( predictions[name][i][-1].start, predictions[name][i][-1].stop ) )
#                         del predictions[name][i][-1]
#                     else:
#                         break
#                 if len( predictions[name][i] ) == 0:
#                     del predictions[name][i]
#                 else:
#                     i += 1
#             logging.debug( "=> %s" % ( str( predictions[name] ) ) )
#             if len( predictions[name] ) == 0:
#                 del predictions[name]


def _create_protplotdata(fname, idx, values):
    """
    create the data for the protein plot
    @param fname the file to write into
    """
    f = open(fname, "w")
    for name in idx:
        tab = []
        for rf in idx[name]:
            for i in list(idx[name][rf].keys()):
                p = idx[name][rf][i]

                tab.append(
                    (
                        values[p].get_position(),
                        values[p].get_score(),
                        values[p].get_query(),
                        values[p].get_yreuq(),
                    )
                )

        # sort by position (and reverse score)
        tab.sort(key=lambda a: (a[0], -a[1]))
        # remove entries at the same position with smaller score
        for i in reversed(list(range(1, len(tab)))):
            if tab[i][0] == tab[i - 1][0] and tab[i][1] <= tab[i - 1][1]:
                del tab[i]

        if len(tab) > 0:
            f.write("%d %s %d %f %f %f\n" % (tab[0][0] - 1, name, rf, 0, 0, 0))
            li = tab[0][0] - 1
            for i in range(len(tab)):
                if li != tab[i][0] - 1:
                    f.write("%d %s %d %f %f %f\n" % (li + 1, name, rf, 0, 0, 0))
                    f.write("%d %s %d %f %f %f\n" % (tab[i][0] - 1, name, rf, 0, 0, 0))

                f.write(
                    "%d %s %d %f %f %f\n"
                    % (tab[i][0], name, rf, tab[i][1], tab[i][2], tab[i][3])
                )
                li = tab[i][0]

            f.write("%d %s %d %f %f %f\n" % (tab[-1][0] + 1, name, rf, 0, 0, 0))
        del tab

    f.close()


def _nonoverlapping_greedy(
    hits, idx, values, finovl, clipfac, scoresel="e", pavg=None, havg=None
):
    """
    greedy determination of a non-overlapping subset of the given hits
    @param[in] hits
    @param scoresel: what should be the basis of the score at a position. 'e':evalue, 'b':bitscore, 'h':height
    @param pavg: score at a position is the average if true, othewise the sum
    @param havg: score of an initial hit is the average score of the position's scores if true, else: sum
    @param finovl: maximum overlap between predictions (nt)
    @param clipfac

    @return the subset
    """

    #    pergen
    #    hitspg = {}
    #    for i in range( len( hits ) ):
    #        try:
    #            hitspg[ hits[i].name ].append( hits[i] )
    #        except KeyError:
    #            hitspg[ hits[i].name ] = [ hits[i] ]
    #    hits = []
    #    for name in hitspg.keys():
    #        hitspg[ name ].sort( key = lambda x: ( x.score, x.stop - x.start ) )
    #
    #    while True:
    #        curbest = []
    #        for name in hitspg.keys():
    #            if len( hitspg[name] ) > 0:
    #                curbest.append( hitspg[ name ][-1] )
    #                del hitspg[ name ][-1]
    #            else:
    #                del hitspg[ name ]
    #        if len( curbest ) == 0:
    #            break
    #        else:
    #            curbest.sort( key = lambda x: ( x.score, x.stop - x.start ), reverse = True )
    #
    #        hits += curbest
    cliplb = 1.0 / float(clipfac)
    clipub = 1.0 * float(clipfac)

    # determine (greedily) a set of non overlapping predictions
    hits.sort(key=lambda x: (x.score, x.stop - x.start), reverse=True)
    i = 0
    while i < len(hits):
        j = 0
        while j < i:
            logging.debug("testing %s %s" % (str(hits[i]), str(hits[j])))
            # compute length of intersection and determine if this makes
            # more than 20% of one of the hits
            cap, cup = hits[i].capcup(hits[j], circular=False, size=0)
            #             logging.error( "%s-%s %d %f %f" % ( hits[i].name, hits[j].name, cap, cap / float( hits[i].length( False, 0 ) ), cap / float( hits[j].length( False, 0 ) ) ) )
            if cap > finovl:
                logging.debug("cap %d", cap)
                #                 overlap = True
                break
            j += 1

        # overlap was detected
        # - if it is the same gene and the heights of the initial predictions
        #   are appropriately similar then the predictions are adapted accordingly
        # - otherwise the new hit (the smaller hit is deleted)
        if j < i:  # hits[i].strand == hits[j].strand and \
            if (
                hits[i].name == hits[j].name
                and cliplb <= (hits[i].score / hits[j].score) <= clipub
            ):
                logging.debug("clip overlapping %s %s" % (str(hits[i]), str(hits[j])))
                hits[i], hits[j] = _clip_predictions(
                    hits[i], hits[j], idx[hits[i].name], values
                )
                logging.debug("                 %s %s" % (str(hits[i]), str(hits[j])))

                if hits[i] is None:
                    del hits[i]
                else:
                    i = i + 1
            else:
                logging.debug("\tdelete %s" % (str(hits[i])))
                del hits[i]
        else:
            i = i + 1

    return hits


def _nonoverlapping_greedy_sorted(hits, finovl):
    """
    old version of greedy determination of non overlapping initial prediction
    subset (here the hits are traversed sorted by starting position and
    only neigbouring pairs are checked, this is faster returns worse results
    as conflict occur more easily)
    @param hits: a list of initial predictions
    @param finovl: maximum overlap between predictions
    @return: a non-overlapping subset of hits
    """

    # remove hits that overlap by more than 20% (of the length of one of the
    # features)
    hits.sort(key=lambda x: x.start)
    i = 0
    while i < len(hits):
        k = i + 1

        while k < len(hits):
            # simpe check if i-1 and i can overlap
            if hits[i].stop < hits[k].start:
                break

            # compute length of intersection and determine if this makesstr
            # more than 20% of one of the hits
            cap, cup = hits[i].capcup(hits[k], circular=False, size=0)
            if cap > finovl:
                #                print "cap", cap, "cup", cup, "len(i)", float( hits[i].length( False, 0 ) ), "len(k)", float( hits[k].length( False, 0 ) )
                # print cap / float( hits[i].length( False, 0 ) ), cap / float(
                # hits[k].length( False, 0 ) )
                scorei = (hits[i].score, hits[i].stop - hits[i].start + 1)
                scorek = (hits[k].score, hits[k].stop - hits[k].start + 1)
                if scorei > scorek:
                    logging.debug("favour %s delete %s" % (str(hits[i]), str(hits[k])))
                    del hits[k]
                    k = k - 1
                elif scorek > scorei:
                    logging.debug("favour %s delete %s" % (str(hits[k]), str(hits[i])))
                    del hits[i]
                    i = i - 1
                    break
                else:
                    logging.warning(
                        "found two equaly good overlapping hits: %s %s"
                        % (hits[i - 1], hits[i])
                    )

            k = k + 1

        i = i + 1
    return hits


def _update_hit(hit, idx, values, circular=False, length=0):
    """
    update hit and position values if a hit got extended or shrinked

    for extended hits also add 0 to idx for the new positions

    @param[in,out] hit the hit to update
    @param[in] idx index from (rf,position) to index in values
    @param[in] values quality values
    @param[in] circular genome circularity
    @param[in] length genome length
    """

    # determine the 1st and last position of the old hit
    start = None
    stop = None

    for p in feature.crange(hit.start, hit.stop + 1, 1, circular, length):
        if p in idx[hit.rf]:
            if start is None or p < start.get_position():
                start = values[idx[hit.rf][p]]
            if stop is None or p > stop.get_position():
                stop = values[idx[hit.rf][p]]
        else:
            idx[hit.rf][p] = len(values)
            values.append(position_values(p))
    conspos = []
    for p in feature.crange(hit.start, hit.stop + 1, 1, circular, length):
        #         if p in not idx[hit.rf]:
        #             idx[hit.rf][ p ] = len( values )
        #             values.append( position_values( p ) )

        # set the query and yreuq values to
        # the value of the 1st or last position
        # of the old hit
        if p <= start.get_position():
            values[-1].set_query(start.get_query())
            values[-1].set_yreuq(start.get_yreuq())
        else:
            values[-1].set_query(stop.get_query())
            values[-1].set_yreuq(stop.get_yreuq())

        conspos.append(values[idx[hit.rf][p]])

    hit.update(conspos)

    if hit.score == 0:
        raise UpdateZeroScore(hit)
