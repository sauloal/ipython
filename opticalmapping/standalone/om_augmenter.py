#!/usr/bin/python

import os
import sys

from om_shared import *


def parse_args(args):
    parser = argparse.ArgumentParser(description="Bionano Genomics MAP parser")
    parser.add_argument( 'infile',                                    help="MAP file"                                               )
    parser.add_argument( '-g'    , '--count'  , action='store_false', help="DO NOT perform global count"                            )
    parser.add_argument( '-c'    , '--conf'   , action='store_false', help="DO NOT perform confidence stats"                        )
    
    args    = parser.parse_args(args=args)

    return args

def main(args):
    valid_fields        = gen_valid_fields(valid_fields_g)

    infile              = args.infile
    DO_GLOBAL_COUNT     = args.count
    DO_CONFIDENCE_STATS = args.conf
    oufile              = infile + ".augmented.tsv"

    if not os.path.exists(infile):
        print "input file %s does not exists" % infile
        sys.exit(1)
        
    if os.path.isdir(infile):
        print "input file %s is a folder" % infile
        sys.exit(1)
    
    
    
    print "saving to %s" % oufile

    data, headers, names, seman, types, indexer, groups, ref_maps_from, query_maps_from, filters_csv = parse_file(infile, valid_fields)

    print "NAMES"  , names
    print "TYPES"  , types
    #print "HEADERS", "\n".join( headers )
    #print "DATA" , data[1]
    #print "INDEX", indexer.keys()[0], indexer[indexer.keys()[0]]
    
    print "file has %5d maps and %3d chromosomes" % (len(indexer["QryContigID"]), len(indexer["RefContigID"]))

    
    if DO_GLOBAL_COUNT:
        print "PRINTING GLOBAL COUNT"
        for RefContigID in sorted(groups["RefContigID_QryContigID"]):
            print "chromosome %2d has %4d maps" % ( RefContigID, len(groups["RefContigID_QryContigID"][RefContigID]))
        print


    
    if DO_CONFIDENCE_STATS:
        print "PRINTING CONFIDENCE STATS"
        for QryContigID in sorted(groups["QryContigID_RefContigID"]):
            print "query %5d maps to %2d chromosomes" % (QryContigID, len(groups["QryContigID_RefContigID"][QryContigID]))
            XmapEntryIDs   = groups["QryContigID_XmapEntryID"][QryContigID].keys()
            Confidences    = [groups["XmapEntryID_Confidence"][x].keys()[0] for x in XmapEntryIDs]
            
            print " confidences         ", Confidences
            max_confidence = max(Confidences)
            
            print " max confidence      ", max_confidence
            print " max confidence chrom", data[list(groups["XmapEntryID_Confidence"][XmapEntryIDs[Confidences.index(max_confidence)]][max_confidence])[0]][seman["RefContigID"]]
        print
    

    
    print "CREATING REPORT:", oufile 
    data      = [ KeyedTuple(x, labels=names)._asdict() for x in data ]
    
    with open(oufile, "w") as reporter:
        reporter.write("\n".join(headers[:-2]) + "\n#\n")
        reporter.write("# FIELDS:\n")
        reporter.write(        "\n".join( [ "# %-39s: %s" % ( x, valid_fields['helps_t'][x] ) for x in valid_fields['names'  ] ] ) + "\n#\n")
        
        reporter.write("#h " + "\t".join( [ "%-39s"       % ( x                             ) for x in valid_fields['names'  ] ] ) + "\n"  )
        reporter.write("#f " + "\t".join( [ "%-39s"       % (    valid_fields['types'  ][x] ) for x in valid_fields['names'  ] ] ) + "\n"  )

        for RefContigID in sorted(groups["RefContigID_QryContigID"]):
            QryContigIDs = groups["RefContigID_QryContigID"][RefContigID]
            
            for QryContigID in sorted(QryContigIDs):
                data_poses           = list(groups["RefContigID_QryContigID"][RefContigID][QryContigID])
                all_data_poses       = list(indexer["QryContigID"][QryContigID])
                data_vals            = [ data[x] for x in data_poses ]
    
                stats                = stats_from_data_vals(RefContigID, QryContigID, groups, indexer, data, data_vals, all_data_poses)
    
                #print "RefContigID %4d QryContigID %6d" % ( RefContigID, QryContigID )
                for data_val in data_vals:
                    cigar                = data_val["HitEnum"]
                    cigar_matches, cigar_insertions, cigar_deletions = process_cigar(cigar)
                    
                    Alignment            = data_val["Alignment"]
                    alignment_count_queries, alignment_count_refs, alignment_count_refs_colapses, alignment_count_queries_colapses = process_alignment(Alignment)
                    
                    for stat in stats:
                        data_val[stat] = stats[stat]
    
                    data_val["_meta_alignment_count_queries"              ] = alignment_count_queries
                    data_val["_meta_alignment_count_queries_colapses"     ] = alignment_count_refs_colapses
                    data_val["_meta_alignment_count_refs"                 ] = alignment_count_refs
                    data_val["_meta_alignment_count_refs_colapses"        ] = alignment_count_queries_colapses
                
                    data_val["_meta_cigar_deletions"                      ] = cigar_deletions
                    data_val["_meta_cigar_insertions"                     ] = cigar_insertions
                    data_val["_meta_cigar_matches"                        ] = cigar_matches

                    data_val["_meta_proportion_query_len_gapped"          ] = (data_val['_meta_len_qry_match_gapped'] * 1.0)/ data_val["QryLen"]
                    data_val["_meta_proportion_query_len_no_gap"          ] = (data_val['_meta_len_qry_match_no_gap'] * 1.0)/ data_val["QryLen"] 
    
                    #print " ", " ".join( ["%s %s" % (x, str(data_val[x])) for x in sorted(data_val)] )
                    reporter.write(            "\t".join( [ str(data_val[x])                                   for x in valid_fields['names'  ] ] ) + "\n"  )



if __name__ == '__main__':
    if len(sys.argv) ==1:
        print "no arguments given"
        sys.exit(1)

    args         = parse_args(sys.argv[1:])

    main(args)
    
"""    
# $ cd D:\Plextor\data\Acquisitie\BioNanoGenomics\MyLycopersicumWorkspace_31022015\Imports; C:\Program Files\BioNano Genomics\RefAligner\WindowsRefAligner.exe -f -ref D:\Plextor\data\Acquisitie\BioNanoGenomics\MyLycopersicumWorkspace_31022015\Imports\S_lycopersicum_chromosomes.2.50.BspQI-BbvCI.cmap -i D:\Plextor\data\Acquisitie\BioNanoGenomics\MyLycopersicumWorkspace_31022015\Imports\EXP_REFINEFINAL1.cmap -o S_lycopersicum_chromosomes.2.50.BspQI-BbvCI_to_EXP_REFINEFINAL1 -endoutlier 1e-2 -outlier 1e-4 -extend 1 -FN 0.08 -FP 0.8 -sf 0.2 -sd 0 -sr 0.02 -res 2.9 -resSD 0.7 -mres 2.0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 2 -T 1e-6 -stdout -stderr

# r3498 $Header: http://svn.bnm.local:81/svn/informatics/RefAligner/branches/3480/RefAligner.cpp 3470 2014-12-17 19:29:21Z tanantharaman $

# FLAGS: USE_SSE=0 USE_AVX=0 USE_MIC=0 USE_PFLOAT=1 USE_RFLOAT=1 DEBUG=1 VERB=1
# XMAP File Version:    0.2
# Label Channels:       1
# Reference Maps From:  S_lycopersicum_chromosomes.2.50.BspQI-BbvCI_to_EXP_REFINEFINAL1_r.cmap
# Query Maps From:      S_lycopersicum_chromosomes.2.50.BspQI-BbvCI_to_EXP_REFINEFINAL1_q.cmap
#h XmapEntryID  QryContigID     RefContigID     QryStartPos     QryEndPos       RefStartPos     RefEndPos       Orientation     Confidence      HitEnum QryLen  RefLen  LabelChannel    Alignment
#f int          int             int             float           float           float           float           string          float           string  float   float   int             string
1       141     1       528400.6        571697.5        10672   54237.5 +       6.65    4M2D2M  1439123.5       21805821        1       "(1,34)(2,34)(3,35)(4,36)(5,37)(6,38)(8,38)(9,39)"
2       174     1       21236.5 1568390 10672   1553561 +       79.35   2M3D1M1D1M1D4M1I2M1D2M1D1M2I2D9M3I3M1D6M1D2M2D1M1D6M1D1M1D1M2D2M2D1M1I1D1M1D5M2D4M2D1M2D2M1D2M1D3M1D1M1D2M3I3D1M1D1M3D2M3D1M2I1D1M2D1M1D1M1I2D3M2I1M1D2M1D1M1D1M2I3D3M3D1M2D1M1D1M1D5M2D12M     1568410 21805821        1       "(1,2)(2,2)(3,3)(6,4)(7,4)(9,5)(11,6)(12,7)(13,8)(14,9)(15,11)(16,12)(18,13)(19,14)(20,15)(21,15)(24,18)(25,19)(26,20)(27,21)(28,22)(29,23)(30,24)(31,25)(32,26)(33,30)(34,31)(35,32)(37,33)(38,34)(39,35)(40,36)(41,37)(42,38)(44,39)(45,40)(47,41)(48,41)(50,42)(51,43)(52,44)(53,45)(54,46)(55,47)(57,48)(59,49)(60,50)(62,50)(63,51)(66,52)(68,54)(69,55)(70,55)(71,56)(72,57)(73,58)(74,59)(76,60)(77,60)(78,61)(79,62)(80,63)(82,64)(83,64)(86,65)(87,66)(89,67)(90,68)(92,69)(93,70)(94,71)(95,72)(96,72)(98,73)(99,74)(103,78)(105,79)(109,80)(110,81)(111,82)(114,82)(116,85)(119,86)(120,87)(121,87)(124,89)(125,90)(126,91)(127,94)(128,95)(129,95)(130,96)(132,97)(134,98)(138,101)(139,102)(140,103)(143,104)(144,104)(146,105)(147,105)(149,106)(151,107)(152,108)(153,109)(154,110)(155,111)(158,112)(159,113)(160,114)(161,115)(162,116)(163,117)(164,118)(165,119)(166,120)(167,121)(168,122)(169,123)"
"""