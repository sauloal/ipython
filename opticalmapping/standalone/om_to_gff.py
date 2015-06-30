#!/usr/bin/python

import os
import sys

from om_shared import *


def parse_args(args):
    parser = argparse.ArgumentParser(description="Bionano Genomics augmented MAP to GFF converter")
    parser.add_argument( 'infile',                                                                    help="AUGMENTED file"                                 )
    parser.add_argument( '-x'    , '--exclude-cols',                                  action='store', help="Exclude column from GFF" )
    parser.add_argument( '-z'    , '--exclude-cols-from-file',                        action='store', help="File containing names of columns to exclude from GFF" )
    parser.add_argument( '-n'    , '--names',                                         action='store', help="Names of reference chromosome. Eg: Ch0|Ch1|Ch3 Ch0,Ch1,Ch2 Ch0:Ch1:Ch2" )
    parser.add_argument( '-f'    , '--names-from-file',                               action='store', help="File containing names of reference chromosome. One per line" )
    parser.add_argument( '-s'    , '--sep'  , '--separator',  default=",",                            help="Separator for chromosome names. Eg: | , :" )
    
    ##genome-build source buildName
    ##species NCBI_Taxonomy_URI
    
    args    = parser.parse_args(args=args)

    if args.names is None and args.names_from_file is None:
        print "either --names or --names-from-file has to be defined"
        sys.exit(1)
    
    if args.names_from_file is not None:
        if not os.path.exists(args.names_from_file):
            print "names file %s does not exists" % args.names_from_file
            sys.exit(1)
            
        if os.path.isdir(args.names_from_file):
            print "names file %s is a folder" % args.names_from_file
            sys.exit(1)

        args.names = []
        with open(args.names_from_file, 'r') as fhd:
            for line in fhd:
                line = line.strip()
                
                if len(line) == 0:
                    continue
                
                if line[0] == "#":
                    continue
                
                args.names.append(line)
                
    elif args.names is not None:
        args.names = args.names.split( args.sep )


    
    if args.exclude_cols_from_file is not None:
        if not os.path.exists(args.exclude_cols_from_file):
            print "columns to exclude file %s does not exists" % args.exclude_cols_from_file
            sys.exit(1)
            
        if os.path.isdir(args.exclude_cols_from_file):
            print "columns to exclude file %s is a folder" % args.exclude_cols_from_file
            sys.exit(1)

        args.exclude_cols = []
        with open(args.exclude_cols_from_file, 'r') as fhd:
            for line in fhd:
                line = line.strip()
                
                if len(line) == 0:
                    continue
                
                if line[0] == "#":
                    continue
                
                args.exclude_cols.append(line)
                
    elif args.exclude_cols is not None:
        args.exclude_cols = args.exclude_cols.split( args.sep )

    else:
        args.exclude_cols = []

    return args

def main(args):
    valid_fields        = gen_valid_fields(valid_fields_g)
    infile              = args.infile
    chromosome_names    = args.names
    exclude_cols        = args.exclude_cols
    oufile              = infile + ".gff3"
    source_name         = "IrysView"
    feature_name_full   = "optical_contig"
    feature_name_piece  = "optical_contig_piece"

    feature_name_full1  = "gene"
    feature_name_full2  = "mRNA"
    feature_name_piece  = "CDS"
    
    id_prefix           = "om_"


    if not os.path.exists(infile):
        print "input file %s does not exists" % infile
        sys.exit(1)
        
    if os.path.isdir(infile):
        print "input file %s is a folder" % infile
        sys.exit(1)
    
    print "saving to %s" % oufile

    data, headers, names, seman, types, indexer, groups, ref_maps_from, query_maps_from, filters_csv = parse_file(infile, valid_fields)

    filters = gen_filter(filters_csv, valid_fields)

    print "NAMES"  , names
    #print "HEADERS", "\n".join( headers )
    print "TYPES"  , types
    #print "DATA" , data[1]
    #print "INDEX", indexer.keys()[0], indexer[indexer.keys()[0]]
    
    print "file has %5d maps and %3d chromosomes" % (len(indexer["QryContigID"]), len(indexer["RefContigID"]))

    assert len(indexer["RefContigID"]       ) <= len(chromosome_names), "number of chromosome differ from %d to %d\n%s\n%s" % (len(indexer["RefContigID"]       ), len(chromosome_names),     indexer["RefContigID"].keys() , chromosome_names)
    assert max(indexer["RefContigID"].keys()) <= len(chromosome_names), "number of chromosome differ from %d to %d"         % (max(indexer["RefContigID"].keys()), len(chromosome_names))

    print chromosome_names
    
    print "CREATING GFF: ", oufile
    with open(oufile, "w") as fhd:
        linecount = 0
        #fhd.write("/home/assembly/nobackup/mummer/MUMmer3.23/1502/solanum_lycopersicum_heinz/SL2.40ch12.fa /home/assembly/nobackup/mummer/MUMmer3.23/1502/solanum_pennellii_scaffold/final.assembly.fasta\n")
        
        fhd.write("##gff-version 3\n")
        
        for RefContigID in sorted(groups["RefContigID_RefStartPos"]):
            RefStartPoses = groups["RefContigID_RefStartPos"][RefContigID]
            RefEndPoses   = groups["RefContigID_RefEndPos"  ][RefContigID]
            ref_min_pos   = min( [min(RefEndPoses), min(RefStartPoses)] )
            ref_max_pos   = max( [max(RefEndPoses), max(RefStartPoses)] )
            fhd.write("##sequence-region %s %d %d\n" % (chromosome_names[RefContigID-1], int(ref_min_pos), int(ref_max_pos)))
        
        
        fhd.write( "#\n" + "\n".join([ "# XMAP "+x[1:] for x in headers[:-2] ]) + "\n#\n")

        
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
                        done_QryContigID[QryContigID][RefContigID] = True
                    
                    qry_rows             = list(groups["RefContigID_QryContigID"][RefContigID][QryContigID])
                    
                    ref_lens             = [ ( data[x]["RefStartPos"], data[x]["RefEndPos"] ) for x in qry_rows ]
                    qry_lens             = [ ( data[x]["QryStartPos"], data[x]["QryEndPos"] ) for x in qry_rows ]
        
                    ref_no_gap_len       = sum( [ max(x)-min(x) for x in ref_lens ] )
                    ref_min_coord        = min( [ min(x)        for x in ref_lens ] )
                    ref_max_coord        = max( [ max(x)        for x in ref_lens ] )
                    
                    qry_no_gap_len       = sum( [ max(x)-min(x) for x in qry_lens ] )
                    qry_min_coord        = min( [ min(x)        for x in qry_lens ] )
                    qry_max_coord        = max( [ max(x)        for x in qry_lens ] )


                    chromosome_name = chromosome_names[RefContigID-1]

                    attributes_keys_G = [
                        [ 'ID'  , "%s%s_%d"   % ( id_prefix, chromosome_name, QryContigID ) ],
                        [ 'Name', "%s%s_%d"   % ( id_prefix, chromosome_name, QryContigID ) ]
                    ]

                    for filter_data in filters:
                        attributes_keys_G.append( [ "_meta_filter_"+filter_data[0].lower(), filter_data[1] + '_' + str(filter_data[3]) ] )
        
                    attributes_G  = ";".join("=".join([k,str(v)]) for k,v in attributes_keys_G)
                    
                    line_G = [ chromosome_name, source_name, feature_name_full1, int(ref_min_coord), int(ref_max_coord), '.', '.', '.',  attributes_G]
                    fhd.write(      "\t".join([str(x) for x in line_G]) + "\n")




                    attributes_keys_G = [
                        [ 'ID'    , "%s%s_%d_m" % ( id_prefix, chromosome_name, QryContigID ) ],
                        [ 'Name'  , "%s%s_%d_m" % ( id_prefix, chromosome_name, QryContigID ) ],
                        [ 'Parent', "%s%s_%d"   % ( id_prefix, chromosome_name, QryContigID ) ]
                    ]

                    for filter_data in filters:
                        attributes_keys_G.append( [ "_meta_filter_"+filter_data[0].lower(), filter_data[1] + '_' + str(filter_data[3]) ] )
        
                    attributes_G  = ";".join("=".join([k,str(v)]) for k,v in attributes_keys_G)
                    
                    line_G = [ chromosome_name, source_name, feature_name_full2, int(ref_min_coord), int(ref_max_coord), '.', '.', '.',  attributes_G]
                    fhd.write(      "\t".join([str(x) for x in line_G]) + "\n")




                    qry_num = 1
                    for qry_row_pos in qry_rows:
                        qry_row         = data[qry_row_pos]
                        
                        RefStartPos     = qry_row["RefStartPos"]
                        RefEndPos       = qry_row["RefEndPos"  ]
                        RefLen          = qry_row["RefLen"     ]
                        
                        QryStartPos     = qry_row["QryStartPos"]
                        QryEndPos       = qry_row["QryEndPos"  ]
                        QryLen          = qry_row["QryLen"     ]
                        
                        Confidence      = qry_row["Confidence" ]
                        Orientation     = qry_row["Orientation"]
                        HitEnum         = qry_row["HitEnum"    ]
            
                        attributes_keys = [
                            [ 'ID'    , "%s%s_%d_m_%06d_c" % ( id_prefix, chromosome_name, QryContigID, qry_num ) ],
                            [ 'Name'  , "%s%s_%d_m_%06d_c" % ( id_prefix, chromosome_name, QryContigID, qry_num ) ],
                            [ 'Parent', "%s%s_%d_m"        % ( id_prefix, chromosome_name, QryContigID          ) ],
                            [ 'Gap'   , HitEnum                                  ],
                        ]
                        
                        for k in sorted(qry_row):
                            if k in exclude_cols:
                                continue
                            attributes_keys.append( [ k.lower(), qry_row[k] ] )
            
                        for filter_data in filters:
                            attributes_keys.append( [ "_meta_filter_"+filter_data[0].lower(), filter_data[1] + '_' + str(filter_data[3]) ] )
            
                        attributes  = ";".join("=".join([k,str(v)]) for k,v in attributes_keys)

                        #http://www.ensembl.org/info/website/upload/gff.html
                        #        seqname          source       feature             start             end             score       strand       frame attribute 
                        line = [ chromosome_name, source_name, feature_name_piece, int(RefStartPos), int(RefEndPos), Confidence, Orientation, '.',  attributes]
                        fhd.write(      "\t".join([str(x) for x in line]) + "\n")

                        qry_num += 1





        #for RefContigID in sorted(groups["RefContigID_RefStartPos"]):
        #    RefStartPoses = groups["RefContigID_RefStartPos"][RefContigID]
        #
        #    for RefStartPosG in sorted(RefStartPoses):
        #        pos_rows  = list(RefStartPoses[RefStartPosG])
        #        
        #        for pos_row_pos in pos_rows:
        #            pos_row         = data[pos_row_pos]
        #
        #            QryContigID     = pos_row["QryContigID"]
        #            
        #            RefStartPos     = pos_row["RefStartPos"]
        #            RefEndPos       = pos_row["RefEndPos"  ]
        #            RefLen          = pos_row["RefLen"     ]
        #            
        #            QryStartPos     = pos_row["QryStartPos"]
        #            QryEndPos       = pos_row["QryEndPos"  ]
        #            QryLen          = pos_row["QryLen"     ]
        #            
        #            Confidence      = pos_row["Confidence" ]
        #            Orientation     = pos_row["Orientation"]
        #            HitEnum         = pos_row["HitEnum"    ]
        #
        #            attributes_keys = [
        #                [ 'ID'  , QryContigID ],
        #                [ 'Name', QryContigID ],
        #                [ 'Gap' , HitEnum     ],
        #            ]
        #            
        #            for k in sorted(pos_row):
        #                if k in exclude_cols:
        #                    continue
        #                attributes_keys.append( [ k.lower(), pos_row[k] ] )
        #
        #            for filter_data in filters:
        #                attributes_keys.append( [ "_meta_filter_"+filter_data[0].lower(), filter_data[1] + '_' + str(filter_data[3]) ] )
        #
        #            attributes  = ";".join("=".join([k,str(v)]) for k,v in attributes_keys)
        #
        #            #http://www.ensembl.org/info/website/upload/gff.html
        #            #        seqname                          source       feature       start             end             score       strand       frame attribute 
        #            line = [ chromosome_names[RefContigID-1], source_name, feature_name, int(RefStartPos), int(RefEndPos), Confidence, Orientation, '.',  attributes]
        #            fhd.write(      "\t".join([str(x) for x in line]) + "\n")

        print



if __name__ == '__main__':
    if len(sys.argv) ==1:
        print "no arguments given"
        sys.exit(1)

    args         = parse_args(sys.argv[1:])

    main(args)