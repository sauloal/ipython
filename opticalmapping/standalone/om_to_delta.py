#!/usr/bin/python

import os
import sys

from om_shared import *


def parse_args(args):
    parser = argparse.ArgumentParser(description="Bionano Genomics augmented MAP to Mummerplot Delta converter")
    parser.add_argument( '-r'    , '--reference',     action='store', help="reference fasfa file" )
    parser.add_argument( '-q'    , '--query'    ,     action='store', help="query fasfa file"     )
    parser.add_argument( 'infile',                                    help="AUGMENTED file"                                         )
    
    args    = parser.parse_args(args=args)

    return args

def main(args):
    valid_fields        = gen_valid_fields(valid_fields_g)
    infile              = args.infile
    oufile              = infile + ".delta"

    if not os.path.exists(infile):
        print "input file %s does not exists" % infile
        sys.exit(1)
        
    if os.path.isdir(infile):
        print "input file %s is a folder" % infile
        sys.exit(1)
    
    print "saving to %s" % oufile

    data, headers, names, seman, types, indexer, groups, ref_maps_from, query_maps_from, filters_csv = parse_file(infile, valid_fields)

    print "NAMES"  , names
    #print "HEADERS", "\n".join( headers )
    print "TYPES"  , types
    #print "DATA" , data[1]
    #print "INDEX", indexer.keys()[0], indexer[indexer.keys()[0]]
    
    print "file has %5d maps and %3d chromosomes" % (len(indexer["QryContigID"]), len(indexer["RefContigID"]))

    

    print "CREATING DELTA: ", oufile
    with open(oufile, "w") as fhd:
        linecount = 0
        #fhd.write("/home/assembly/nobackup/mummer/MUMmer3.23/1502/solanum_lycopersicum_heinz/SL2.40ch12.fa /home/assembly/nobackup/mummer/MUMmer3.23/1502/solanum_pennellii_scaffold/final.assembly.fasta\n")
        if args.reference is not None and args.query is not None:
            fhd.write("%s %s\n" % (args.reference, args.query))
        else:
            fhd.write("%s %s\n" % (ref_maps_from, query_maps_from))
        
        fhd.write("NUCMER\n")

        sum_ref_len = 0
        sum_qry_len = 0
        
        data = [ KeyedTuple(x, labels=names)._asdict() for x in data ]
        
        done_QryContigID = {}
        for RefContigID in sorted(groups["RefContigID_RefStartPos"]):
            RefStartPoses = groups["RefContigID_RefStartPos"][RefContigID]
            RefLen        = 0

            for RefStartPosG in sorted(RefStartPoses):
                pos_rows = list(RefStartPoses[RefStartPosG])
                
                for pos_row_pos in pos_rows:
                    pos_row = data[pos_row_pos]

                    QryContigID = pos_row["QryContigID"]
                    
                    if QryContigID not in done_QryContigID:
                        done_QryContigID[QryContigID] = {}
                    
                    first_time = True
                    
                    if RefContigID in done_QryContigID[QryContigID]:
                        continue
                    
                    else:
                        if len(done_QryContigID[QryContigID]) == 0:
                            done_QryContigID[QryContigID][RefContigID] = sum_qry_len

                        else:
                            done_QryContigID[QryContigID][RefContigID] = done_QryContigID[QryContigID][done_QryContigID[QryContigID].keys()[0]]
                            first_time = False
                    
                    sum_qry_len_local = done_QryContigID[QryContigID][RefContigID]
                    
                    line_count        = 0
                    qry_rows          = list(groups["RefContigID_QryContigID"][RefContigID][QryContigID])
                    
                    for qry_row_pos in qry_rows:
                        qry_row     = data[qry_row_pos]
                        
                        RefStartPos = qry_row["RefStartPos"] + sum_ref_len
                        RefEndPos   = qry_row["RefEndPos"  ] + sum_ref_len
                        RefLen      = qry_row["RefLen"     ]
                        
                        QryStartPos = qry_row["QryStartPos"] + sum_qry_len_local
                        QryEndPos   = qry_row["QryEndPos"  ] + sum_qry_len_local
                        QryLen      = qry_row["QryLen"     ]
                        
                        num_errors  = 0
                        sim_erros   = 0
                        stop_codons = 0
                        
                        header      = [int(RefContigID), int(QryContigID), int(RefLen     ), int(QryLen   ) ]
                        line        = [int(RefStartPos), int(RefEndPos  ), int(QryStartPos), int(QryEndPos), num_errors, sim_erros, stop_codons]

                        if line_count == 0:
                            fhd.write(">" + " ".join([str(x) for x in header]) + "\n")

                        fhd.write(      " ".join([str(x) for x in line  ]) + "\n0\n")

                        line_count += 1
                    
                    if first_time:
                        sum_qry_len += QryLen
            
            sum_ref_len += RefLen
            #break
        print



if __name__ == '__main__':
    if len(sys.argv) ==1:
        print "no arguments given"
        sys.exit(1)

    args         = parse_args(sys.argv[1:])

    main(args)