#!/usr/bin/python

import os
import sys

from om_shared import *


def parse_args(args):
    parser = argparse.ArgumentParser(description="Bionano Genomics augmented MAP filter")
    parser.add_argument( 'infile',                                    help="AUGMENTED file"                                             )
    parser.add_argument( '-l'    , '--list'   , action='store_true' , help="List Fields and Operators"                                  )
    parser.add_argument( '-f'    , '--filter' , action='append'     , help="Filters [Field:Function(%s):Value]" % ", ".join(sorted(valid_operators.keys())))
    
    args    = parser.parse_args(args=args)

    return args


def main(args):
    valid_fields = gen_valid_fields(valid_fields_g)
    infile       = args.infile

    if args.list:
        print "LIST OF FIELDS"
        print "", "\n ".join( ["%-41s: %-6s : %s"% (valid_field_name, valid_fields['types'  ][valid_field_name], valid_fields['helps'  ][valid_field_name]) for valid_field_name in valid_fields['names'  ]] )
        print "LIST OF OPERATORS"
        print "", "\n ".join( [ "%-9s: %s" % (x,valid_operators[x][1]) for x in sorted(valid_operators) ] )
        sys.exit(0)

    if not os.path.exists(infile):
        print "input file %s does not exists" % infile
        sys.exit(1)
        
    if os.path.isdir(infile):
        print "input file %s is a folder" % infile
        sys.exit(1)
    
    
    filters = gen_filter(args.filter, valid_fields)
    
    oufile = infile
    for field_name, field_operator_name, field_operator, field_value_str, field_value in filters:
        oufile += ('_' + '_'.join( [ field_name, field_operator_name, field_value_str ] )).replace(',', '_')
    
    print "saving to %s" % oufile

    data, headers, names, seman, types, indexer, groups, ref_maps_from, query_maps_from, filters_csv = parse_file(infile, valid_fields)
    data = [ KeyedTuple(x, labels=names)._asdict() for x in data ]

    print "NAMES"  , names
    #print "HEADERS", "\n".join( headers )
    print "TYPES"  , types
    #print "DATA" , data[1]
    #print "INDEX", indexer.keys()[0], indexer[indexer.keys()[0]]
    
    print "file has %5d maps and %3d chromosomes" % (len(indexer["QryContigID"]), len(indexer["RefContigID"]))

    
    print "CREATING REPORT:", oufile + ".report.tsv"
    with open(oufile + ".report.tsv", "w") as reporter:
        reporter.write("\n".join(headers[:-2]) + "\n#\n")
        
        reporter.write("# FILTERS:\n")
        for field_name, field_operator_name, field_operator, field_value_str, field_value in filters:
            reporter.write( "# FILTER : %-39s: %3s : %s\n" % ( field_name, field_operator_name, field_value_str ) )
        reporter.write( "\n\n" )
        
        reporter.write("#h " + "\t".join( [ "%-39s"       % ( x                             ) for x in valid_fields['names'  ] ] ) + "\n")
        reporter.write("#f " + "\t".join( [ "%-39s"       % (    valid_fields['types'  ][x] ) for x in valid_fields['names'  ] ] ) + "\n")



        for RefContigID in sorted(groups["RefContigID_QryContigID"]):
            QryContigIDs = groups["RefContigID_QryContigID"][RefContigID]
            
            for QryContigID in sorted(QryContigIDs):
                data_poses            = list(groups["RefContigID_QryContigID"][RefContigID][QryContigID])
                data_vals             = []
                valid_data_poses      = []
                
                all_data_poses = list(indexer["QryContigID"][QryContigID])
                for data_pos in list(all_data_poses):
                    data_val   = data[data_pos]
                    filter_res = all([ field_operator( data_val[field_name], field_value ) for field_name, field_operator_name, field_operator, field_value_str, field_value in filters])
                    if filter_res:
                        if data_pos in data_poses:
                            data_vals.append(data_val)
                        valid_data_poses.append(data_pos)
    
                if len(data_vals) > 0:
                    stats = stats_from_data_vals(RefContigID, QryContigID, groups, indexer, data, data_vals, valid_data_poses)
        
                    #print "RefContigID %4d QryContigID %6d" % ( RefContigID, QryContigID )
                    for data_val in data_vals:
                        for stat in stats:
                            data_val[stat] = stats[stat]

                        filter_res = all([ field_operator( data_val[field_name], field_value ) for field_name, field_operator_name, field_operator, field_value_str, field_value in filters])
                        if not filter_res:
                            continue

                        #print " ", " ".join( ["%s %s" % (x, str(data_val[x])) for x in sorted(data_val)] )
                        reporter.write(            "\t".join( [ str(data_val[x])                                   for x in valid_fields['names'  ] ] ) + "\n")
    print




if __name__ == '__main__':
    if len(sys.argv) ==1:
        print "no arguments given"
        sys.exit(1)

    args         = parse_args(sys.argv[1:])

    main(args)